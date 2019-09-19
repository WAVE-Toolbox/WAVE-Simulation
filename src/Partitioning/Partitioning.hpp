#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <cmath>
#include <vector>

#ifdef USE_GEOGRAPHER
#include <geographer/ParcoRepart.h>
#include <geographer/AuxiliaryFunctions.h>
#ifdef USE_GEOGRAPHER_WRAPPERS
#include <geographer/Wrappers.h>
#endif
#endif
namespace KITGPI
{
    //! \brief Partitioning namespace
    namespace Partitioning
    {

        /*! \brief inter node distribution define the grid topology by sizes NX, NY, and NZ from configuration   
               *    Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1
            \param config configuration object
            \param commShot communicator of a shot domain
            */
        template <typename ValueType>
        dmemo::DistributionPtr gridPartition(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commShot)
        {
            common::Grid3D grid(config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<IndexType>("NX"));
            // distribute the grid onto available processors
            return (std::make_shared<dmemo::GridDistribution>(grid, commShot));
        }

#ifdef USE_GEOGRAPHER
        /*! \brief 
            \param config configuration object
            \param commShot communicator of a shot domain
            \param coords std::vector of the coordinate Densevectors (x,y,z)
            \param graph CSR Matrix of the graph (combined derivative matrices)
            \param weights DenseVector of node weights for each gridpoint
            */
        template <typename ValueType>
        dmemo::DistributionPtr graphPartition(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commShot, std::vector<scai::lama::DenseVector<ValueType>> &coords, scai::lama::CSRSparseMatrix<ValueType> &graph, scai::lama::DenseVector<ValueType> &weights, ITI::Tool tool=ITI::Tool::geoKmeans)
        {
            std::string dimension = config.get<std::string>("dimension");

            IndexType dimensions = 0;
            // transform to lower cases
            std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
            // Assert correctness of input values
            SCAI_ASSERT_ERROR(dimension.compare("3d") == 0 || dimension.compare("2d") == 0, "Unkown dimension");
            if (dimension.compare("2d") == 0) {
                dimensions = 2;
                coords.pop_back();
            } else {
                dimensions = 3;
            }

            struct ITI::Settings settings;
            settings.dimensions = dimensions;
            settings.noRefinement = true;
            settings.verbose = false;
            settings.minBorderNodes = 100;
            settings.multiLevelRounds = 9;
            settings.numBlocks = commShot->getSize();
            settings.coarseningStepsBetweenRefinement = 1;

			//const IndexType N = graph.getNumRows();
			const IndexType M = graph.getNumValues(); 
			//set minimum gain to 1% of the average number of edges per block
			settings.minGainForNextRound = M/settings.numBlocks*0.01;    
            //settings.maxKMeansIterations = 10;

            //settings.minSamplingNodes = -1;
            settings.storeInfo = true;
		    //look at ITI Settings.h for the enum ITI::Tool
            settings.initialPartition = tool;
           
            if( tool==ITI::Tool::geographer ){
            	settings.noRefinement = false;
            	settings.initialPartition = ITI::Tool::geoKmeans;
            }
            
            //if hierarchical kmeans, read hierarchy levels
			if( tool==ITI::Tool::geoHierKM or tool==ITI::Tool::geoHierRepart ){
				std::string hL= config.get<std::string>("hierLevels");
				std::vector<std::string> vhl = ITI::aux<IndexType,ValueType>::split( hL, ',');
				for( unsigned int i=0; i<vhl.size(); i++)
					settings.hierLevels.push_back( std::stoi(vhl[i]) );
			}

			//if multisection, read cuts per dimension
			if( tool==ITI::Tool::geoMS ){
				std::string cpd = config.get<std::string>("cutsPerDim");
				std::vector<std::string> vcpd = ITI::aux<IndexType,ValueType>::split( cpd, ',');
				for( unsigned int i=0; i<vcpd.size(); i++){
					settings.cutsPerDim.push_back( std::stoi(vcpd[i]) );
				}

				SCAI_ASSERT_EQ_ERROR( settings.cutsPerDim.size(), (unsigned long) settings.dimensions, "Dimensions and cuts per dimensions must agree" );
			}

            //in case we do local refinement, we change the edge weights to 1
            //if( not settings.noRefinement)
            //update, 30/04: change the weights anyway to get a correct cut
            {
            	//change all edge weights to 1
				scai::lama::CSRStorage<ValueType>& localStorage = graph.getLocalStorage();
				scai::hmemo::HArray<ValueType> localValues = localStorage.getValues();

				for( unsigned int i=0; i<localValues.size(); i++ ){
					localValues[i] = 1.0;
				}
				
				scai::hmemo::HArray<IndexType> localIA = localStorage.getIA();
				scai::hmemo::HArray<IndexType> localJA = localStorage.getJA();
				localStorage.swap( localIA, localJA, localValues);
			}

			try{
				settings.mappingRenumbering = config.get<bool>("mappingRenumbering");
			}catch(...){
				settings.mappingRenumbering = false;
			}

            ITI::Metrics<ValueType> metrics(settings); //by default, settings.numBlocks = p (where p is: mpirun -np p ...)

            if (commShot->getRank() == 0) {
                settings.print(std::cout);
            }

            std::vector<lama::DenseVector<ValueType>> weightVector;
            weightVector.push_back(weights);

			scai::lama::DenseVector<IndexType> partition;

			if( ITI::to_string(tool).rfind("geo",0)==0 ){
				partition = ITI::ParcoRepart<IndexType, ValueType>::partitionGraph(graph, coords, weightVector, settings, metrics);
			}else{
#ifdef USE_GEOGRAPHER_WRAPPERS
				bool nodeWeightsUse = true; //usign unit weights
		    	partition = ITI::Wrappers<IndexType,ValueType>::partition( graph, coords, weightVector, nodeWeightsUse, tool, settings, metrics );
#endif	    	
    		}

            if (config.get<bool>("partitionWrite"))
                partition.writeToFile(config.get<std::string>("partitionFilename") + ".mtx");

            dmemo::DistributionPtr dist = scai::dmemo::generalDistributionByNewOwners(partition.getDistribution(), partition.getLocalValues());

            //redistribute all data to get metrics
            scai::dmemo::DistributionPtr noDistPtr(new scai::dmemo::NoDistribution(graph.getNumRows()));
            graph.redistribute(dist, noDistPtr);
            partition.redistribute(dist);
            weightVector[0].redistribute(dist);

            metrics.getAllMetrics(graph, partition, weightVector, settings);

            if (commShot->getRank() == 0) {
                metrics.print(std::cout);
            }

            return (dist);
        }
#endif
        /*! \brief calculation of the weights for the absorbing boundary
            \param config configuration object
            \param dist distributionPtr of the model
            \param modelCoordinates coordinate object
            */
        template <typename ValueType>
        scai::lama::DenseVector<ValueType> BoundaryWeights(Configuration::Configuration const &config, dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType weight)
        {

            std::string dimension = config.get<std::string>("dimension");
            std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

            auto BoundaryWidth = config.get<IndexType>("BoundaryWidth");

            hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
            dist->getOwnedIndexes(ownedIndexes);

            lama::VectorAssembly<ValueType> assembly;
            assembly.reserve(ownedIndexes.size());
            if (config.get<IndexType>("DampingBoundary")) {
                for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
                    Acquisition::coordinate3D coordinatedist = modelCoordinates.edgeDistance(coordinate);

                    scai::IndexType min = 0;
                    if (coordinatedist.x < coordinatedist.y) {
                        min = coordinatedist.x;
                    } else {
                        min = coordinatedist.y;
                    }

                    if (dimension.compare("3d") == 0) {
                        min = coordinatedist.min();
                    }

                    if (config.get<IndexType>("FreeSurface") == 0) {
                        if (min < BoundaryWidth) {
                            assembly.push(ownedIndex, weight);
                        }
                    } else {
                        IndexType HorizontalMin = 0;
                        if (dimension.compare("3d") == 0) {
                            HorizontalMin = !((coordinatedist.x) < (coordinatedist.z)) ? (coordinatedist.z) : (coordinatedist.x);
                        } else {
                            HorizontalMin = coordinatedist.x;
                        }

                        if (coordinate.y < BoundaryWidth) {
                            if (HorizontalMin < BoundaryWidth) {
                                assembly.push(ownedIndex, weight);
                            }
                        } else if (min < BoundaryWidth) {
                            assembly.push(ownedIndex, weight);
                        }
                    }
                }
            }

            lama::DenseVector<ValueType> weights;
            weights.allocate(dist);
            weights = 1.0;
            weights.fillFromAssembly(assembly);
            // weights.setContextPtr(ctx);
            if (config.get<bool>("weightsWrite"))
                weights.writeToFile(config.get<std::string>("weightsFilename") + ".mtx");

            return (weights);
        }
    }
}
