#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/dmemo/GridDistribution.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration/Configuration.hpp"

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"
#include "Acquisition/suHandler.hpp"

#include "ForwardSolver/ForwardSolver.hpp"

#include "ForwardSolver/Derivatives/DerivativesFactory.hpp"
#include "ForwardSolver/ForwardSolverFactory.hpp"
#include "Modelparameter/ModelparameterFactory.hpp"
#include "Wavefields/WavefieldsFactory.hpp"

#include "Common/HostPrint.hpp"

#include "CheckParameter/CheckParameter.hpp"

#include <ParcoRepart.h>


using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA

    common::Settings::parseArgs(argc, argv);

   // typedef float ValueType;
    double start_t, end_t; /* For timing */

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    // following lines should be in a extra function in check parameter class
    IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    IndexType firstSnapshot = 0, lastSnapshot = 0, incSnapshot = 0;
    if (config.get<IndexType>("snapType") > 0) {
        firstSnapshot = static_cast<IndexType>(config.get<ValueType>("tFirstSnapshot") / config.get<ValueType>("DT") + 0.5);
        lastSnapshot = static_cast<IndexType>(config.get<ValueType>("tLastSnapshot") / config.get<ValueType>("DT") + 0.5);
        incSnapshot = static_cast<IndexType>(config.get<ValueType>("tIncSnapshot") / config.get<ValueType>("DT") + 0.5);
    }

    IndexType partitionedOut = static_cast<IndexType>(config.get<ValueType>("PartitionedOut"));

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    IndexType npS = config.get<IndexType>("ProcNS");
    IndexType npX = config.get<IndexType>("ProcNX");
    IndexType npY = config.get<IndexType>("ProcNY");
    IndexType npZ = config.get<IndexType>("ProcNZ");

    // following lines should be part of checkParameter.tpp
//     if (commAll->getSize() != npS * npX * npY * npZ) {
//         HOST_PRINT(commAll, "\n Error: Number of MPI processes (" << commAll->getSize() << ") doesn't match the number of processes specified in " << argv[1] << ": ProcNS * ProcNX * ProcNY * ProcNZ = " << npS * npX * npY * npZ << "\n")
//         return (2);
//     }

    //number of processes inside a shot domain
    IndexType npNpS=commAll->getSize()/npS;
    
    // Build subsets of processors for the shots

    common::Grid2D procAllGrid(npS, npNpS);
    IndexType procAllGridRank[2];

    procAllGrid.gridPos(procAllGridRank, commAll->getRank());

    // communicator for set of processors that solve one shot

    dmemo::CommunicatorPtr commShot = commAll->split(procAllGridRank[0]);

    // this communicator is used for reducing the solutions of problems

    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    // inter node distribution
    // define the grid topology by sizes NX, NY, and NZ from configuration
    // Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1

    common::Grid3D grid(config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<IndexType>("NX"));

    common::Grid3D procGrid(npY, npZ, npX);
    // distribute the grid onto available processors
    
   //  dmemo::DistributionPtr dist(new dmemo::GridDistribution(grid, commShot, procGrid));
   int size=config.get<IndexType>("NY")*config.get<IndexType>("NZ")*config.get<IndexType>("NX");
   dmemo::DistributionPtr dist(new dmemo::BlockDistribution(size,commShot));

  
    // Create an object of the mapping (3D-1D) class Coordinates

    Acquisition::Coordinates<ValueType> modelCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"), dist, ctx);

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::VectorAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size()); 
    
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        Acquisition::coordinate3D coordinatedist = modelCoordinates.edgeDistance(coordinate);

        scai::IndexType min = 0;
        if (coordinatedist.x < coordinatedist.y) {
            min = coordinatedist.x;
        } else {
            min = coordinatedist.y;
        }
        
        if (min < config.get<IndexType>("BoundaryWidth")) {
  //      if (coordinatedist.min() <  config.get<IndexType>("BoundaryWidth")){
        assembly.push(ownedIndex, 1.5);
        }
    }
    lama::DenseVector<ValueType> weights;
    weights.allocate(dist);
    weights=1.0;
    weights.fillFromAssembly(assembly);
    weights.setContextPtr(ctx);

    weights.writeToFile("weights.mtx");
    
    verbose = config.get<bool>("verbose");
    HOST_PRINT(commAll, "\nSOFI" << dimension << " " << equationType << " - LAMA Version\n\n");

    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    derivatives->init(dist, ctx, config, modelCoordinates, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    
    
    /* --------------------------------------- */
    /* Call partioner                          */
    /* --------------------------------------- */
    IndexType dimensions=0;
   // transform to lower cases
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("3d") == 0 || dimension.compare("2d") == 0, "Unkown dimension");

    if (dimension.compare("2d") == 0) {
        dimensions=2;
    }
    if (dimension.compare("3d") == 0) {
        dimensions=3;
    }

    HOST_PRINT(commAll, modelCoordinates.getCoordinates()[2] << "\n\n");
    HOST_PRINT(commAll, weights << "\n\n");
    HOST_PRINT(commAll, derivatives->getDxf() << "\n\n");
    scai::lama::CSRSparseMatrix<ValueType> graph=derivatives->getCombinedMatrix();
    HOST_PRINT(commAll, graph << "\n\n");
    HOST_PRINT(commAll, dimensions << "\n\n");
    auto coords=modelCoordinates.getCoordinates();
    
    if (dimensions==2){
    coords.pop_back();
    }
    
//std::cout << graph.getNumValues() << std::endl;

if( commShot->getRank() ==0 ){
	std::cout << graph.getLocalStorage().getValues()[0] << std::endl;
	//std::cout << graph.getLocalStorage().getValues().sum() << std::endl;
	std::cout << graph.getLocalStorage().getValues()[100] << std::endl;
}

//coords[0].writeToFile("coordsX.mtx");
//coords[1].writeToFile("coordsY.mtx");

//change all edge weights to 1
/*
{
	CSRStorage<ValueType>& localStorage = graph.getLocalStorage();
	//scai::hmemo::WriteAccess<ValueType> values(localStorage.getValues());

	scai::hmemo::HArray<IndexType> localIA = localStorage.getIA();
	scai::hmemo::HArray<IndexType> localJA = localStorage.getJA();
	scai::hmemo::HArray<ValueType> localValues = localStorage.getValues();

	for( unsigned int i=0; i<localValues.size(); i++ ){
		localValues[i] = 1;
	}

	localStorage.swap( localIA, localJA, localValues);
}
*/

PRINT(commShot->getRank() << ": " << (DenseVector<ValueType>(graph.getLocalStorage().getValues())).sum() );	
std::string graphFile = "gpi100x100_edge_weights.graph";
std::string coordsFile = "gpi100x100.graph.xyz";
if( commShot->getRank()==0 ) PRINT( "writing graph and coords to files " << graphFile << " and " << coordsFile );

//ITI::FileIO<IndexType, ValueType>::writeCoords( coords, coordsFile );
ITI::FileIO<IndexType, ValueType>::writeGraph( graph, graphFile, 1 /*for edges weights*/ );

    struct Settings settings;
    settings.dimensions=dimensions;
    settings.noRefinement = false;
    settings.verbose = false;
    settings.minBorderNodes = 100;
    settings.multiLevelRounds = 6;
   	settings.numBlocks=commShot->getSize();
    settings.coarseningStepsBetweenRefinement = 3;
    //settings.maxKMeansIterations = 10;
    //settings.minSamplingNodes = -1;
    settings.writeInFile = false;
    //settings.bisect = true;
    settings.initialPartition = InitialPartitioningMethods::KMeans;

    struct Metrics metrics(settings);      //by default, settings.numBlocks = p (where p is: mpirun -np p ...)
    
    if( commShot->getRank() ==0 ){
		settings.print( std::cout );
	}

    scai::lama::DenseVector<IndexType> partition = ITI::ParcoRepart<IndexType,ValueType>::partitionGraph(graph, coords, weights, settings, metrics);
	    
    partition.writeToFile("partitition.mtx");
    
    dmemo::DistributionPtr distFromPartition = scai::dmemo::generalDistributionByNewOwners( partition.getDistribution(), partition.getLocalValues() );
    //    dmemo::DistributionPtr distFromPartition(new dmemo::GridDistribution(grid, commShot, procGrid));
     //  dmemo::DistributionPtr distFromPartition(new dmemo::BlockDistribution(size,commShot));
       
    
    derivatives->redistributeMatrices(distFromPartition);
    
   	//redistribute all data to get metrics
    scai::dmemo::DistributionPtr noDistPtr( new scai::dmemo::NoDistribution( graph.getNumRows() ));
 //TODO: if we remove the partition.redistribute works but I am not sure the the computation of the cut is correct:
 // local refinement says it gives a gain but the final cut is more than the first cut
 // if we remove any other, the computeCut hangs
   	graph.redistribute( distFromPartition , noDistPtr );
    partition.redistribute( distFromPartition );
    weights.redistribute( distFromPartition );

	metrics.getAllMetrics(graph, partition, weights, settings);

	if( commShot->getRank() ==0 ){
		metrics.print( std::cout );
	}

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, distFromPartition);

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Sources<ValueType> sources(config, modelCoordinates, ctx, distFromPartition);
    CheckParameter::checkSources<ValueType>(config, sources, commShot);
    Acquisition::Receivers<ValueType> receivers;
    if (!config.get<bool>("useReceiversPerShot")) {
        receivers.init(config, modelCoordinates, ctx, distFromPartition);
        CheckParameter::checkReceivers<ValueType>(config, receivers, commShot);
    }

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, distFromPartition);
    model->prepareForModelling(modelCoordinates, ctx, distFromPartition, commShot);
    CheckParameter::checkNumericalArtefeactsAndInstabilities<ValueType>(config, *model, commShot);

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */

    HOST_PRINT(commAll, "", "ForwardSolver ...\n")
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));
    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
    HOST_PRINT(commAll, "", "ForwardSolver prepared\n")

    dmemo::BlockDistribution shotDist(sources.getNumShots(), commInterShot);

    for (IndexType shotNumber = shotDist.lb(); shotNumber < shotDist.ub(); shotNumber++) {
        /* Update Source */
        if (!config.get<bool>("runSimultaneousShots"))
            sources.init(config, modelCoordinates, ctx, distFromPartition, shotNumber);
        if (config.get<bool>("useReceiversPerShot")) {
            receivers.init(config, modelCoordinates, ctx, distFromPartition, shotNumber);
            CheckParameter::checkReceivers<ValueType>(config, receivers, commShot);
        }

        HOST_PRINT(commShot, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n",
                   "Total Number of time steps: " << tStepEnd << "\n");
        wavefields->resetWavefields();

        start_t = common::Walltime::get();

        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

            if (tStep % 100 == 0 && tStep != 0) {
                HOST_PRINT(commShot, "Calculating time step " << tStep << "\n");
            }

            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);

            if (config.get<IndexType>("snapType") > 0 && tStep >= firstSnapshot && tStep <= lastSnapshot && (tStep - firstSnapshot) % incSnapshot == 0) {
                wavefields->write(config.get<IndexType>("snapType"), config.get<std::string>("WavefieldFileName"), tStep, *derivatives, *model, partitionedOut);
            }
        }

        end_t = common::Walltime::get();
        HOST_PRINT(commShot, "Finished time stepping (shot " << shotNumber << ") in " << end_t - start_t << " sec.\n");

        receivers.getSeismogramHandler().normalize();

        if (!config.get<bool>("runSimultaneousShots")) {
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
        } else {
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename"), modelCoordinates);
        }

        solver->resetCPML();
    }
    return 0;
}
