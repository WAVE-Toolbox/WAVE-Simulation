#pragma once

#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../IO/IO.hpp"

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <cmath>
#include <vector>

#include <scai/dmemo/RedistributePlan.hpp>

#ifdef USE_GEOGRAPHER
#include <geographer/AuxiliaryFunctions.h>
#include <geographer/ParcoRepart.h>
#ifdef USE_GEOGRAPHER_WRAPPERS
#include <geographer/Wrappers.h>
#endif
#endif

#include "../Common/Common.hpp"
#include <scai/common/ContextType.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/partitioning/Partitioning.hpp>

namespace KITGPI
{
    //! \brief Partitioning namespace
    namespace Partitioning
    {

        //         /*! \brief inter node distribution define the grid topology by sizes NX, NY, and NZ from configuration
        //             *Attention: LAMA uses row-major indexing while WAVE-Simulation-3D uses column-major, so switch dimensions, x-dimension has stride 1 z-dimension has stride 2 y-dimension has stride 3
        //             \param config configuration object
        //             \param commShot communicator of a shot domain
        //             */
        // template <typename ValueType>
        IndexType getShotDomain(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commAll)
        {

            /* Definition of shot domains */
            IndexType shotDomainDefinition = config.get<int>("ShotDomainDefinition");
            IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
            IndexType shotDomain = 0;                                           // will contain the domain to which this processor belongs
            if (shotDomainDefinition == 0) {
                // Definition by number of shot domains

                /* Definition of shot domains */

                IndexType npDomain = commAll->getSize() / numShotDomains; // number of processors for each shot domain

                if (commAll->getSize() != numShotDomains * npDomain) {
                    HOST_PRINT(commAll, "\n Error: Number of MPI processes (" << commAll->getSize()
                                                                              << ") is not multiple of shot domains in the configuration: numShotDomains = " << numShotDomains << "\n\n\n")
                    COMMON_THROWEXCEPTION("Number of MPI processes is not multiple of shot domains in the configuration: numShotDomains");
                }

                shotDomain = commAll->getRank() / npDomain;
            } else if (shotDomainDefinition == 1) {
                // All processors on one node build one domain

                shotDomain = commAll->getNodeId();
            } else {
                bool set = common::Settings::getEnvironment(shotDomain, "DOMAIN");

                if (!set) {
                    std::cout << *commAll << ", node = " << commAll->getNodeName()
                              << ", node rank = " << commAll->getNodeRank() << " of " << commAll->getNodeSize()
                              << ": environment variable DOMAIN not set" << std::endl;
                }

                set = commAll->all(set); // make sure that all processors will terminate

                if (!set) {
                    COMMON_THROWEXCEPTION("set is false for process with Rank " << commAll->getRank());
                }
            }
            return (shotDomain);
        }

        /*! \brief inter node distribution define the grid topology by sizes NX, NY, and NZ from configuration   
               *    Attention: LAMA uses row-major indexing while WAVE-Simulation-3D uses column-major, so switch dimensions, x-dimension has stride, 1 z-dimension has stride 2 y-dimension has stride 3
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
        dmemo::DistributionPtr graphPartition(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commShot, std::vector<scai::lama::DenseVector<ValueType>> &coords, scai::lama::CSRSparseMatrix<ValueType> &graph, scai::lama::DenseVector<ValueType> &weights, ITI::Tool tool = ITI::Tool::geoKmeans)
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
            settings.minGainForNextRound = M / settings.numBlocks * 0.01;
            //settings.maxKMeansIterations = 10;

            //settings.minSamplingNodes = -1;
            settings.storeInfo = true;
            //look at ITI Settings.h for the enum ITI::Tool
            settings.initialPartition = tool;

            if (tool == ITI::Tool::geographer) {
                settings.noRefinement = false;
                settings.initialPartition = ITI::Tool::geoKmeans;
            }

            //if hierarchical kmeans, read hierarchy levels
            if (tool == ITI::Tool::geoHierKM or tool == ITI::Tool::geoHierRepart) {
                std::string hL = config.get<std::string>("hierLevels");
                std::vector<std::string> vhl = ITI::aux<IndexType, ValueType>::split(hL, ',');
                for (unsigned int i = 0; i < vhl.size(); i++)
                    settings.hierLevels.push_back(std::stoi(vhl[i]));
            }

            //if multisection, read cuts per dimension
            if (tool == ITI::Tool::geoMS) {
                std::string cpd = config.get<std::string>("cutsPerDim");
                std::vector<std::string> vcpd = ITI::aux<IndexType, ValueType>::split(cpd, ',');
                for (unsigned int i = 0; i < vcpd.size(); i++) {
                    settings.cutsPerDim.push_back(std::stoi(vcpd[i]));
                }

                SCAI_ASSERT_EQ_ERROR(settings.cutsPerDim.size(), (unsigned long)settings.dimensions, "Dimensions and cuts per dimensions must agree");
            }

            //in case we do local refinement, we change the edge weights to 1
            //if( not settings.noRefinement)
            //update, 30/04: change the weights anyway to get a correct cut
            {
                //change all edge weights to 1
                scai::lama::CSRStorage<ValueType> &localStorage = graph.getLocalStorage();
                scai::hmemo::HArray<ValueType> localValues = localStorage.getValues();

                for (scai::IndexType i = 0; i < localValues.size(); i++) {
                    localValues[i] = 1.0;
                }

                scai::hmemo::HArray<IndexType> localIA = localStorage.getIA();
                scai::hmemo::HArray<IndexType> localJA = localStorage.getJA();
                localStorage.swap(localIA, localJA, localValues);
            }

            try {
                settings.mappingRenumbering = config.get<bool>("mappingRenumbering");
            } catch (...) {
                settings.mappingRenumbering = false;
            }

            ITI::Metrics<ValueType> metrics(settings); //by default, settings.numBlocks = p (where p is: mpirun -np p ...)

            if (commShot->getRank() == 0) {
                settings.print(std::cout);
            }

            std::vector<lama::DenseVector<ValueType>> weightVector;
            weightVector.push_back(weights);

            scai::lama::DenseVector<IndexType> partition;

            if (ITI::to_string(tool).rfind("geo", 0) == 0) {
                partition = ITI::ParcoRepart<IndexType, ValueType>::partitionGraph(graph, coords, weightVector,commShot, settings, metrics);
            } else {
#ifdef USE_GEOGRAPHER_WRAPPERS
                bool nodeWeightsUse = true; //usign unit weights
                partition = ITI::Wrappers<IndexType, ValueType>::partition(graph, coords, weightVector, nodeWeightsUse, tool, settings, metrics);
#endif
            }

            dmemo::DistributionPtr dist = scai::dmemo::generalDistributionByNewOwners(partition.getDistribution(), partition.getLocalValues());

            //redistribute all data to get metrics (uncommend for debugging or monitiring the partitioner results)
            //             scai::dmemo::DistributionPtr noDistPtr(new scai::dmemo::NoDistribution(graph.getNumRows()));
            //             graph.redistribute(dist, noDistPtr);
            //             partition.redistribute(dist);
            //             weightVector[0].redistribute(dist);
            //
            //             metrics.getEasyMetrics(graph, partition, weightVector, settings);
            //
            //             if (commShot->getRank() == 0) {
            //                 metrics.print(std::cout);
            //             }

            return (dist);
        }
#endif

        /*! \brief calculation of the node weights (variable fd order + pml)
            \param config configuration object
            \param dist distributionPtr of the model
            \param modelCoordinates coordinate object
            */
        template <typename ValueType>
        scai::lama::DenseVector<ValueType> Weights(Configuration::Configuration const &config, dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {
            SCAI_REGION("KITGPI.Weights")
            // Weights of single operations; valid for all cases
            ValueType MatrixVector2ndOrderWeight = 1.00;
            ValueType VectorAssignmentWeight = 0.25;
            ValueType VectorPlusVectorWeight = 0.43;
            ValueType PMLWeight = 1.5;

            IndexType NumMatrixVector = 0;
            IndexType NumVectorAssignement = 0;
            IndexType NumVectorPlusVector = 0;
            IndexType NumPMLPerDim = 0;

            std::string dimension = config.get<std::string>("dimension");
            std::string type = config.get<std::string>("equationType");

            // transform to lower cases
            std::transform(type.begin(), type.end(), type.begin(), ::tolower);
            std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

            // Assert correctness of input values
            SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
            SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0 || type.compare("sh") == 0, "Unkown type");

            //Number of operations during the time stepping
            // 2D
            if (dimension.compare("2d") == 0 && type.compare("acoustic") == 0) {
                NumMatrixVector = 4;
                NumVectorAssignement = 7;
                NumVectorPlusVector = 0;
                NumPMLPerDim = 2;
            }
            if (dimension.compare("2d") == 0 && type.compare("elastic") == 0) {
                NumMatrixVector = 8;
                NumVectorAssignement = 20;
                NumVectorPlusVector = 0;
                NumPMLPerDim = 4;
            }
            if (dimension.compare("2d") == 0 && type.compare("visco") == 0) {
                NumMatrixVector = 8;
                NumVectorAssignement = 41;
                NumVectorPlusVector = 8;
                NumPMLPerDim = 4;
            }
            if (dimension.compare("2d") == 0 && type.compare("sh") == 0) {
                NumMatrixVector = 4;
                NumVectorAssignement = 7;
                NumVectorPlusVector = 0;
                NumPMLPerDim = 2;
            }

            // 3D
            if (dimension.compare("3d") == 0 && type.compare("acoustic") == 0) {
                NumMatrixVector = 6;
                NumVectorAssignement = 10;
                NumVectorPlusVector = 0;
                NumPMLPerDim = 2;
            }
            if (dimension.compare("3d") == 0 && type.compare("elastic") == 0) {
                NumMatrixVector = 18;
                NumVectorAssignement = 34;
                NumVectorPlusVector = 3;
                NumPMLPerDim = 6;
            }
            if (dimension.compare("3d") == 0 && type.compare("visco") == 0) {
                NumMatrixVector = 18;
                NumVectorAssignement = 72;
                NumVectorPlusVector = 22;
                NumPMLPerDim = 6;
            }

            //runtime of Vector Operations are not influenced by the FDorder;
            ValueType constantWeight = NumVectorAssignement * VectorAssignmentWeight + NumVectorPlusVector * VectorPlusVectorWeight;
            ValueType referenceTotalWeight = NumMatrixVector * MatrixVector2ndOrderWeight + constantWeight;

            hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
            dist->getOwnedIndexes(ownedIndexes);

            lama::DenseVector<ValueType> fdWeights;
            fdWeights.setSameValue(dist, referenceTotalWeight);

            bool useNodeWeights = 1;
            try {
                useNodeWeights = config.get<bool>("useNodeWeights");
            } catch (...) {
                //do nothing... use default useNodeWeights=1. useNodeWeights is only used for debugging
            }

            if (useNodeWeights) {
                lama::VectorAssembly<ValueType> assembly;
                assembly.reserve(ownedIndexes.size());

                std::vector<scai::IndexType> spatialFDorderVec;
                if (config.get<bool>("useVariableFDoperators")) {
                    unsigned int column = 2;
                    Common::readColumnFromFile(config.get<std::string>("gridConfigurationFilename"), spatialFDorderVec, column);
                }

                // Weights of single Matrix vector Products with FDorder 2-12
                ValueType FDWeights[6] = {MatrixVector2ndOrderWeight, 1.56, 2.06, 2.54, 3.00, 3.70};

                //loop over all (local) indeces
                for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

                    const auto layer = modelCoordinates.getLayer(coordinate);
                    IndexType FDOrder = 0;
                    if (config.get<bool>("useVariableFDoperators")) {
                        FDOrder = spatialFDorderVec[layer];
                    } else {
                        FDOrder = config.get<IndexType>("spatialFDorder");
                    }

                    //FDOrder Weights

                    ValueType fdWeight = (NumMatrixVector * FDWeights[FDOrder / 2 - 1] + constantWeight);

                    assembly.push(ownedIndex, fdWeight);
                }

                fdWeights.fillFromAssembly(assembly);
            }

            lama::DenseVector<ValueType> pmlWeights;
            pmlWeights.setSameValue(dist, 0.0);

            if ((config.get<IndexType>("DampingBoundary") == 2) && (useNodeWeights)) {
                lama::VectorAssembly<ValueType> assembly;
                assembly.reserve(ownedIndexes.size());
                auto BoundaryWidth = config.get<IndexType>("BoundaryWidth");
                for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
                    Acquisition::coordinate3D coordinatedist = modelCoordinates.edgeDistance(coordinate);

                    const auto layer = modelCoordinates.getLayer(coordinate);

                    ValueType pmlWeight = 0;

                    IndexType width = std::ceil((float)BoundaryWidth / modelCoordinates.getDHFactor(layer));
                    IndexType xDist = coordinatedist.x / modelCoordinates.getDHFactor(layer);
                    IndexType yDist = coordinatedist.y / modelCoordinates.getDHFactor(layer);
                    IndexType zDist = coordinatedist.z / modelCoordinates.getDHFactor(layer);

                    if (xDist < width) {
                        pmlWeight += NumPMLPerDim * PMLWeight;
                    }

                    if (yDist < width) {
                        IndexType yCoord = coordinate.y / modelCoordinates.getDHFactor(layer);
                        if (yCoord < width) {
                            if (config.get<IndexType>("FreeSurface") == 0) {
                                pmlWeight += NumPMLPerDim * PMLWeight;
                            }
                        } else {
                            pmlWeight += NumPMLPerDim * PMLWeight;
                        }
                    }
                    if (dimension.compare("3d") == 0) {
                        if (zDist < width) {
                            pmlWeight += NumPMLPerDim * PMLWeight;
                        }
                    }

                    assembly.push(ownedIndex, pmlWeight);
                }

                pmlWeights.fillFromAssembly(assembly);
            }

            lama::DenseVector<ValueType> weights = lama::eval<lama::DenseVector<ValueType>>(fdWeights + pmlWeights);
            weights /= referenceTotalWeight;
            weights /= weights.sum();

            //             if (config.get<bool>("weightsWrite")) {
            //                 IO::writeVector(weights, config.get<std::string>("weightsFilename"), config.get<IndexType>("fileFormat"));
            //             }
            return (weights);
        }

        template <typename ValueType>
        dmemo::DistributionPtr graphPartition(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::CommunicatorPtr commShot, scai::dmemo::DistributionPtr BlockDist, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {
            SCAI_REGION("KITGPI.graphPartitionAll")

#ifdef USE_GEOGRAPHER
            double start_t, end_t; /* For timing */
            start_t = common::Walltime::get();
            //             typedef common::shared_ptr<Context> ContextPtr;
            //
            //             ContextPtr loc = Context::getContextPtr( common::context::Host );
            hmemo::ContextPtr loc = hmemo::Context::getContextPtr(scai::common::ContextType::Host);

            HOST_PRINT(commShot, "", "creating partitioner input... \n");

            HOST_PRINT(commShot, "", "caclulate graph for partitioner \n");
            auto &&graph = derivatives.getGraph(BlockDist, modelCoordinates);
            graph.setContextPtr(loc);

            HOST_PRINT(commShot, "", "get coordinate vectors for partitioner \n");
            auto &&coords = modelCoordinates.getCoordinates(BlockDist, ctx);
            coords[0].setContextPtr(loc);
            coords[1].setContextPtr(loc);
            coords[2].setContextPtr(loc);

            HOST_PRINT(commShot, "", "calculate node weights for partitioner \n");
            auto &&weights = Weights(config, BlockDist, modelCoordinates);
            weights.setContextPtr(loc);

            end_t = common::Walltime::get();
            HOST_PRINT(commShot, "", "created partitioner input  in " << end_t - start_t << " sec.\n\n");

            //which tool to use to partition with,

            std::string toolStr;
            try {
                toolStr = config.get<std::string>("graphPartitionTool");
            } catch (...) {
                toolStr = "geoKMeans";
            }

            //get the enum from string
            ITI::Tool tool = ITI::to_tool(toolStr);

            auto dist = KITGPI::Partitioning::graphPartition(config, commShot, coords, graph, weights, tool);

            // derivatives.redistributeMatrices(dist);
            return (dist);
#else
            HOST_PRINT(commShot, "partitioning=2 or useVariableGrid was set, but geographer was not compiled. \nUse < make prog GEOGRAPHER_ROOT= > to compile the partitioner")
            HOST_PRINT(commShot, "\nBlock Distribution will be used instead!\n\n")
            return (BlockDist);
#endif
        }

        template <typename ValueType>
        dmemo::DistributionPtr metisPartition(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::CommunicatorPtr commShot, scai::dmemo::DistributionPtr BlockDist, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {
            SCAI_REGION("KITGPI.metisPartitionAll")

            auto loc = hmemo::Context::getHostPtr();

            HOST_PRINT(commShot, "", "caclulate graph for ParMetis \n");
            auto &&graph = derivatives.getGraph(BlockDist, modelCoordinates);

            HOST_PRINT(commShot, "", "caclulate weights for ParMetis \n");

            auto &&weights = Weights(config, BlockDist, modelCoordinates);

            HOST_PRINT(commShot, "", "create partitioning by ParMetis \n");

            graph.setContextPtr(loc);

            scai::dmemo::DistributionPtr newDist = BlockDist;

            if (scai::partitioning::Partitioning::canCreate("PARMETIS")) {

                auto myPartitioning = scai::partitioning::Partitioning::create("PARMETIS");

                scai::hmemo::HArray<IndexType> newLocalOwners;

                float procWeight = 1.0f;

                scai::hmemo::HArray<float> vertexWeights; // currently only float is supported

                weights.buildLocalValues(vertexWeights);

                myPartitioning->squarePartitioningW(newLocalOwners, graph, vertexWeights, procWeight);

                auto plan = dmemo::redistributePlanByNewOwners(newLocalOwners, BlockDist);

                newDist = plan.getTargetDistributionPtr();
            } else {
                HOST_PRINT(commShot, "", "ATTENTION: PARMETIS partitioning not supported, will use BLOCK distribution")
            }

            return newDist;
        }
    }
}
