#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <cmath>
#include <vector>

#ifdef USE_GEOGRAPHER
#include <ParcoRepart.h>
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
            common::Grid3D procGrid(config.get<IndexType>("ProcNY"), config.get<IndexType>("ProcNZ"), config.get<IndexType>("ProcNX"));
            // distribute the grid onto available processors
            return (std::make_shared<dmemo::GridDistribution>(grid, commShot, procGrid));
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
        dmemo::DistributionPtr graphPartition(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commShot, std::vector<scai::lama::DenseVector<ValueType>> &coords, scai::lama::CSRSparseMatrix<ValueType> &graph, scai::lama::DenseVector<ValueType> &weights)
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

            struct Settings settings;
            settings.dimensions = dimensions;
            settings.noRefinement = true;
            settings.verbose = false;
            settings.minBorderNodes = 100;
            settings.multiLevelRounds = 9;
            settings.numBlocks = commShot->getSize();
            settings.coarseningStepsBetweenRefinement = 1;
            //settings.maxKMeansIterations = 10;
            //settings.minSamplingNodes = -1;
            settings.writeInFile = true;
            settings.initialPartition = InitialPartitioningMethods::KMeans;

            struct Metrics metrics(settings); //by default, settings.numBlocks = p (where p is: mpirun -np p ...)

            if (commShot->getRank() == 0) {
                settings.print(std::cout);
            }

            scai::lama::DenseVector<IndexType> partition = ITI::ParcoRepart<IndexType, ValueType>::partitionGraph(graph, coords, weights, settings, metrics);
            partition.writeToFile("partitition.mtx");
            dmemo::DistributionPtr dist = scai::dmemo::generalDistributionByNewOwners(partition.getDistribution(), partition.getLocalValues());

            //redistribute all data to get metrics
            scai::dmemo::DistributionPtr noDistPtr(new scai::dmemo::NoDistribution(graph.getNumRows()));
            graph.redistribute(dist, noDistPtr);
            partition.redistribute(dist);
            weights.redistribute(dist);

            metrics.getAllMetrics(graph, partition, weights, settings);

            if (commShot->getRank() == 0) {
                metrics.print(std::cout);
            }

            return (dist);
        }
#endif
        /*! \brief calculatio of the weights for the absorbing boundary
            \param config configuration object
            \param dist distributionPtr of the model
            \param modelCoordinates coordinate object
            */
        template <typename ValueType>
        scai::lama::DenseVector<ValueType> BoundaryWeights(Configuration::Configuration const &config, dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates,ValueType weight)
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
            //    weights.setContextPtr(ctx);

            weights.writeToFile("weights.mtx");

            return (weights);
        }
    }
}
