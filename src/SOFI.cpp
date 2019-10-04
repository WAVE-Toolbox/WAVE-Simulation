#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"
#include "Configuration/Configuration.hpp"
#include "Configuration/ValueType.hpp"
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

extern bool verbose; // global variable definition

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA --SCAI_SETTINGS=domains.txt

    common::Settings::parseArgs(argc, argv);

    double start_t, end_t; /* For timing */
    double globalStart_t, globalEnd_t; /* For timing */
    globalStart_t = common::Walltime::get();
    
    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);
    verbose = config.get<bool>("verbose");

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    HOST_PRINT(commAll, "\n SOFI++ " << dimension << " " << equationType << " - LAMA Version\n");
    HOST_PRINT(commAll, "","  - Running on " << commAll->getSize() << " mpi processes -\n\n");
    
    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    std::string settingsFilename; // filename for processor specific settings
    if (common::Settings::getEnvironment(settingsFilename, "SCAI_SETTINGS")) {
        // each processor reads line of settings file that matches its node name and node rank
        common::Settings::readSettingsFile(settingsFilename.c_str(), commAll->getNodeName(), commAll->getNodeRank());
    }

    /* --------------------------------------- */
    /* coordinate mapping (3D<->1D)            */
    /* --------------------------------------- */

    Acquisition::Coordinates<ValueType> modelCoordinates(config);

    if (config.get<bool>("useVariableGrid")) {
        CheckParameter::checkVariableGrid(config, commAll, modelCoordinates);
        for (int layer=0;layer<modelCoordinates.getNumLayers();layer++){
            HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinates.getNGridpoints(layer)); 
        }
        auto numGridpointsRegular=config.get<IndexType>("NX")*config.get<IndexType>("NY")*config.get<IndexType>("NZ");
        HOST_PRINT(commAll, "\n Number of gripoints total: " << modelCoordinates.getNGridpoints());
        HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: "  << (float) modelCoordinates.getNGridpoints()/numGridpointsRegular * 100 << "% \n\n");
    }

    /* --------------------------------------- */
    /* context and communicator for shot parallelisation   */
    /* --------------------------------------- */

    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    IndexType shotDomain = Partitioning::getShotDomain(config, commAll); // will contain the domain to which this processor belongs

    // Build subsets of processors for the shots

    dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);

    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
    SCAI_DMEMO_TASK(commShot)

    /* --------------------------------------- */
    /* Distribution                */
    /* --------------------------------------- */

    dmemo::DistributionPtr dist = nullptr;
    if ((config.get<IndexType>("partitioning") == 0) || (config.get<IndexType>("partitioning") == 2)) {
        //Block distribution = starting distribution for graph partitioner
        dist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
    } else if (config.get<IndexType>("partitioning") == 1) {
        SCAI_ASSERT(!config.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
        dist = Partitioning::gridPartition<ValueType>(config, commShot);
    } else {
        COMMON_THROWEXCEPTION("unknown partioning method");
    }

    if (config.get<bool>("coordinateWrite"))
        modelCoordinates.writeCoordinates(dist, ctx, config.get<std::string>("coordinateFilename"),config.get<IndexType>("FileFormat"));

    /* --------------------------------------- */
    /* Factories                               */
    /* --------------------------------------- */

    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));

    /* --------------------------------------- */
    /* Memory estimation                       */
    /* --------------------------------------- */

    HOST_PRINT(commAll, " ============== Memory Estimation: ===============\n\n")

    ValueType memDerivatives = derivatives->estimateMemory(config, dist, modelCoordinates);
    ValueType memWavefileds = wavefields->estimateMemory(dist);
    ValueType memModel = model->estimateMemory(dist);
    ValueType memSolver = solver->estimateMemory(config, dist, modelCoordinates);
    ValueType memTotal = memDerivatives + memWavefileds + memModel + memSolver;

    HOST_PRINT(commAll, " -  Derivative Matrices \t" << memDerivatives << " MB\n");
    HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefileds << " MB\n");
    HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModel << " MB\n");
    HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolver << " MB\n");
    HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotal << " / " << memTotal / dist->getNumPartitions() << " MB ");
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
    if (numShotDomains > 1)
        HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomains << " shot Domains ): \n " << memTotal * numShotDomains << " MB  ");

    HOST_PRINT(commAll, "\n\n ========================================================================\n\n")

    /* --------------------------------------- */
    /* Call partioner */
    /* --------------------------------------- */
    if (config.get<IndexType>("partitioning") == 2) {
             start_t = common::Walltime::get();
        dist = Partitioning::graphPartition(config, ctx, commShot, dist, *derivatives,modelCoordinates);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished graph partitioning in " << end_t - start_t << " sec.\n\n");
    }

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();

    derivatives->init(dist, ctx, modelCoordinates, commShot);

    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "\n", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    //snapshot of the memory count (freed memory doesn't reduce maxAllocatedBytes())
    // std::cout << "+derivatives "  << hmemo::Context::getHostPtr()->getMemoryPtr()->maxAllocatedBytes() << std::endl;
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;
    Acquisition::readAllSettings<ValueType>(sourceSettings, config.get<std::string>("SourceFilename") + ".txt");
    // build Settings for SU?
    //settings = su.getSourceSettings(shotNumber); // currently not working, expecting a sourceSettings struct and not a vector of sourceSettings structs
    //         su.buildAcqMatrixSource(config.get<std::string>("SourceSignalFilename"), modelCoordinates.getDH());
    //         allSettings = su.getSourceSettingsVec();
    

    Acquisition::Sources<ValueType> sources;

    Acquisition::Receivers<ValueType> receivers;

    if (!config.get<bool>("useReceiversPerShot")) {
        receivers.init(config, modelCoordinates, ctx, dist);
    }
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing Acquisition in " << end_t - start_t << " sec.\n\n");
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    model->init(config, ctx, dist, modelCoordinates);
    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing model in " << end_t - start_t << " sec.\n\n");
    
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    wavefields->init(ctx, dist);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing wavefield in " << end_t - start_t << " sec.\n\n");
    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */


    start_t = common::Walltime::get();
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, config.get<ValueType>("DT"));
    solver->prepareForModelling(*model, config.get<ValueType>("DT"));
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing forward solver in " << end_t - start_t << " sec.\n\n");

    ValueType DT = config.get<ValueType>("DT");
    IndexType tStepEnd = Common::time2index(config.get<ValueType>("T"), DT);

    // calculate vector with unique shot numbers and get number of shots
    std::vector<scai::IndexType> uniqueShotNos;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    CheckParameter::checkSources(sourceSettings,modelCoordinates,commAll);
    IndexType numshots = uniqueShotNos.size();

    /* general block distribution of shot domains accorting to their weights */
    IndexType firstShot = 0;
    IndexType lastShot = numshots - 1;

    float processorWeight = 1.0f;
    common::Settings::getEnvironment(processorWeight, "WEIGHT");
    float domainWeight = commShot->sum(processorWeight);

    if (commShot->getRank() == 0) {
        // master processors of shot domains determine the load distribution
        auto shotDist = dmemo::genBlockDistributionByWeight(numshots, domainWeight, commInterShot);
        firstShot = shotDist->lb();
        lastShot = shotDist->ub();
    }
    commShot->bcast(&firstShot, 1, 0);
    commShot->bcast(&lastShot, 1, 0);

    double end_tInit= common::Walltime::get();
    double tInit=end_tInit-globalStart_t;
    HOST_PRINT(commAll, "", "Finished initializing! in " << tInit << " sec.\n\n");
    
    /* --------------------------------------- */
    /* Loop over shots                        */
    /* --------------------------------------- */

    for (IndexType shotInd = firstShot; shotInd < lastShot; shotInd++) {
        IndexType shotNumber = uniqueShotNos[shotInd];
        /* Update Source */
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
        Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);

        sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

        CheckParameter::checkNumericalArtefeactsAndInstabilities<ValueType>(config, sourceSettingsShot, *model,modelCoordinates,shotNumber);

        bool writeSource_bool;
        try { writeSource_bool = config.get<bool>("writeSource");
        }
        catch (...) {
            writeSource_bool = false;
        }
        if (writeSource_bool) {
            lama::DenseMatrix<ValueType> sourcesignal_out = sources.getsourcesignal();
            KITGPI::IO::writeMatrix(sourcesignal_out, config.get<std::string>("writeSourceFilename") + "_shot_" + std::to_string(shotNumber), config.get<IndexType>("fileFormat"));
        }

        if (config.get<bool>("useReceiversPerShot")) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
        }

        HOST_PRINT(commShot, "Start time stepping for shot " << shotInd << " (shot no: " << shotNumber << "), shotDomain = " << shotDomain << "\n",
                   "\nTotal Number of time steps: " << tStepEnd << "\n");
        wavefields->resetWavefields();

        start_t = common::Walltime::get();

        double start_t2=0.0, end_t2=0.0;

        
        
        /* --------------------------------------- */
        /* Loop over time steps                        */
        /* --------------------------------------- */
        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

            if ((tStep-1)% 100 == 0) {
                 start_t2= common::Walltime::get();
                //HOST_PRINT(commShot, " ", "Calculating time step " << tStep << " in shot  " << shotNumber << "\n");
            }
            
            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);

            if (tStep % 100 == 0 && tStep != 0) {
                 end_t2= common::Walltime::get();
                HOST_PRINT(commShot, "", "Calculated " << tStep << " time steps" << " in shot  " << shotNumber << " at t = " << end_t2 - globalStart_t << "\nLast 100 timesteps calculated in " << end_t2 - start_t2 << " sec. - Estimated runtime (Simulation/total): " << (int) ((tStepEnd/100) * (end_t2 - start_t2)) << " / " << (int) ((tStepEnd/100) * (end_t2 - start_t2) + tInit) << " sec.\n\n");
            }
            
            
            if (config.get<IndexType>("snapType") > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), DT) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT)) % Common::time2index(config.get<ValueType>("tincSnapshot"), DT) == 0) {
                wavefields->write(config.get<IndexType>("snapType"), config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".", tStep, *derivatives, *model, config.get<IndexType>("FileFormat"));
            }
        }

        end_t = common::Walltime::get();
        HOST_PRINT(commShot, "Finished time stepping for shot no: " << shotNumber << " in " << end_t - start_t << " sec.\n","");

        receivers.getSeismogramHandler().normalize();

        receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
        solver->resetCPML();
    }
    globalEnd_t = common::Walltime::get();
    HOST_PRINT(commAll, "\nTotal runtime of SOFI: " << globalEnd_t - globalStart_t << " sec.\nSOFI finished!\n\n");
    return 0;
}
