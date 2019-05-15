#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"
#include "Acquisition/suHandler.hpp"
#include "Configuration/Configuration.hpp"
#include "ForwardSolver/ForwardSolver.hpp"

#include "ForwardSolver/Derivatives/DerivativesFactory.hpp"
#include "ForwardSolver/ForwardSolverFactory.hpp"
#include "Modelparameter/ModelparameterFactory.hpp"
#include "Wavefields/WavefieldsFactory.hpp"

#include "CheckParameter/CheckParameter.hpp"
#include "Common/HostPrint.hpp"
#include "Partitioning/Partitioning.hpp"

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA

    common::Settings::parseArgs(argc, argv);

    typedef double ValueType;
    double start_t, end_t; /* For timing */

    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    Configuration::Configuration config(argv[1]);
    verbose = config.get<bool>("verbose");

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    HOST_PRINT(commAll, "\nSOFI" << dimension << " " << equationType << " - LAMA Version\n\n");
    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    /* --------------------------------------- */
    /* coordinate mapping (3D<->1D)            */
    /* --------------------------------------- */

    Acquisition::Coordinates<ValueType> modelCoordinates(config);

    /* --------------------------------------- */
    /* communicator for shot parallelisation   */
    /* --------------------------------------- */

    IndexType npS = config.get<IndexType>("ProcNS");
    //number of processes inside a shot domain
    IndexType npNpS = commAll->getSize() / npS;

    CheckParameter::checkNumberOfProcesses(config, commAll);

    // Build subsets of processors for the shots
    common::Grid2D procAllGrid(npS, npNpS);
    IndexType procAllGridRank[2];
    procAllGrid.gridPos(procAllGridRank, commAll->getRank());

    // communicator for set of processors that solve one shot
    dmemo::CommunicatorPtr commShot = commAll->split(procAllGridRank[0]);

    // this communicator is used for reducing the solutions of problems
    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */

    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    dmemo::DistributionPtr dist = nullptr;
    dmemo::DistributionPtr blockDist = nullptr;
    if ((config.get<bool>("useVariableGrid")) || (config.get<bool>("useGraphPartitioning"))) {
        //temporaray distribution (starting distribution for graph partitioner)
        blockDist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
    } else {
        dist = Partitioning::gridPartition<ValueType>(config, commShot);
    }

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    if ((config.get<bool>("useVariableGrid")) || (config.get<bool>("useGraphPartitioning"))) {
        derivatives->init(blockDist, ctx, config, modelCoordinates, commShot);
    } else {
        derivatives->init(dist, ctx, config, modelCoordinates, commShot);
    }
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Call partioner (only for variable Grids)*/
    /* --------------------------------------- */
    if ((config.get<bool>("useVariableGrid")) || (config.get<bool>("useGraphPartitioning"))) {

        auto graph = derivatives->getCombinedMatrix();
        auto weights = Partitioning::BoundaryWeights(config, blockDist, modelCoordinates,config.get<ValueType>("BoundaryWeights"));
        auto coords = modelCoordinates.getCoordinates(blockDist, ctx);
        modelCoordinates.writeCoordinates(blockDist, ctx, config.get<std::string>("coordinateFilename"));

        std::string toolStr = config.get<std::string>("graphPartitionTool");
		/*
        
		std::map<std::string,ITI::Tool> toolName = { 
			{"geographer",ITI::Tool::geographer}, {"geoKmeans",ITI::Tool::geoKmeans}, {"geoSFC",ITI::Tool::geoSFC}, {"geoMS", ITI::Tool::geoMS},
			{"geoHierKM", ITI::Tool::geoHierKM}, {"geoHierRepart", ITI::Tool::geoHierRepart},
			{"parMetisSFC",ITI::Tool::parMetisSFC}, {"parMetisGeom",ITI::Tool::parMetisGeom}, {"parMetisGraph",ITI::Tool::parMetisGraph},
			{"zoltanRcb",ITI::Tool::zoltanRCB}, {"zoltanRib",ITI::Tool::zoltanRIB}, {"zoltanMJ",ITI::Tool::zoltanMJ}, {"zoltanHsfc",ITI::Tool::zoltanSFC}
		};

		ITI::Tool tool = toolName[toolStr];
		*/
		//convert string to enum
		//TODO: pass it as string and convert it later to Tool
		ITI::Tool tool = ITI::toTool(toolStr);
		
		start_t = common::Walltime::get();
#ifdef USE_GEOGRAPHER
        dist = Partitioning::graphPartition(config, commShot, coords, graph, weights, tool);
#else
        HOST_PRINT(commAll, "useGraphPartitioning or useVariableGrid was set, but geographer was not compiled. \n Use < make prog GEOGRAPHER_ROOT= > to compile the partitioner\n", "\n")
        return (2);
#endif
		end_t = common::Walltime::get();
        //HOST_PRINT(commShot, "Finished time stepping (shot " << shotNumber << ") in " << end_t - start_t << " sec.\n");
		HOST_PRINT(commShot, "Partitioning time " << end_t - start_t << std::endl) ;

        derivatives->redistributeMatrices(dist);
    }

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Sources<ValueType> sources(config, modelCoordinates, ctx, dist);
    Acquisition::Receivers<ValueType> receivers;
    if (!config.get<bool>("useReceiversPerShot")) {
        receivers.init(config, modelCoordinates, ctx, dist);
    }

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, dist);
    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */

    HOST_PRINT(commAll, "", "ForwardSolver ...\n")
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));
    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
    HOST_PRINT(commAll, "", "ForwardSolver prepared\n")

    ValueType DT = config.get<ValueType>("DT");
    IndexType tStepEnd = Common::time2index(config.get<ValueType>("T"), DT);

    dmemo::BlockDistribution shotDist(sources.getNumShots(), commInterShot);

    for (IndexType shotNumber = shotDist.lb(); shotNumber < shotDist.ub(); shotNumber++) {
        /* Update Source */
        if (!config.get<bool>("runSimultaneousShots"))
            sources.init(config, modelCoordinates, ctx, dist, shotNumber);
        if (config.get<bool>("useReceiversPerShot")) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
        }

        HOST_PRINT(commShot, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\nTotal Number of time steps: " << tStepEnd << "\n");
        wavefields->resetWavefields();

        start_t = common::Walltime::get();

        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

            if (tStep % 100 == 0 && tStep != 0) {
                HOST_PRINT(commShot, "Calculating time step " << tStep << "\n");
            }

            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);

            if (config.get<IndexType>("snapType") > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), DT) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT)) % Common::time2index(config.get<ValueType>("tincSnapshot"), DT) == 0) {
                wavefields->write(config.get<IndexType>("snapType"), config.get<std::string>("WavefieldFileName"), tStep, *derivatives, *model, config.get<IndexType>("PartitionedOut"));
            }
        }

        end_t = common::Walltime::get();
        HOST_PRINT(commShot, "Finished time stepping (shot " << shotNumber << ") in " << end_t - start_t << " sec.\n");
		HOST_PRINT(commShot, "time " << end_t - start_t << std::endl) ;

        receivers.getSeismogramHandler().normalize();

        if (!config.get<bool>("runSimultaneousShots")) {
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
        } else {
            receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename"), modelCoordinates);
        }

        solver->resetCPML();
    }
    std::exit(0); //needed in supermuc
    return 0;
}
