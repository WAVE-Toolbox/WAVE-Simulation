#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
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

#include "config.hpp"

using namespace scai;
using namespace KITGPI;

extern bool verbose; // global variable definition

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA --SCAI_SETTINGS=domains.txt

    common::Settings::parseArgs(argc, argv);

//    typedef float ValueType;
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
    verbose = config.get<bool>("verbose");

    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");

    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    HOST_PRINT(commAll, "\nSOFI" << dimension << " " << equationType << " - LAMA Version\n\n");
    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    std::string settingsFilename;    // filename for processor specific settings

    if ( common::Settings::getEnvironment(settingsFilename, "SCAI_SETTINGS") )
    {
        // each processor reads line of settings file that matches its node name and node rank
        common::Settings::readSettingsFile( settingsFilename.c_str(), commAll->getNodeName(), commAll->getNodeRank() );
    }

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */

    /* execution context */

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    /* Definition of shot domains */

    int shotDomain = config.get<int>("ShotDomain");

    int domain;   // will contain the domain to which this processor belongs

    if (shotDomain == 0)
    {
        // Definition by number of shot domains

       IndexType numDomains = config.get<IndexType>("ProcNS");    // total number of shot domains
       IndexType npDomain = commAll->getSize() / numDomains;      // number of processors for each shot domain

       if (commAll->getSize() != numDomains * npDomain) {
           HOST_PRINT(commAll, "\n Error: Number of MPI processes (" << commAll->getSize()
                               << ") is not multiple of shot domains in " << argv[1] << ": ProcNS = " << numDomains << "\n")
           return (2);
       }
 
       domain = commAll->getRank() / npDomain;
    }
    else if (shotDomain == 1 )
    {
        // All processors on one node build one domain

        domain = commAll->getNodeId();
    }
    else 
    {
        bool set = common::Settings::getEnvironment( domain, "DOMAIN" );

        if (!set)
        {
            std::cout << *commAll << ", node = " << commAll->getNodeName() 
                      << ", node rank = " << commAll->getNodeRank() << " of " << commAll->getNodeSize()
                      << ": environment variable DOMAIN not set" << std::endl;
        }
 
        set = commAll->all( set );   // make sure that all processors will terminate

        if (!set)
        {
            return(2);
        }
    }

    CheckParameter::checkNumberOfProcesses(config, commAll);

    // Build subsets of processors for the shots

    dmemo::CommunicatorPtr commShot = commAll->split(domain);

    // this communicator is used for reducing the solutions of problems

    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());

    SCAI_DMEMO_TASK(commShot)

    // inter node distribution
    // define the grid topology by sizes NX, NY, and NZ from configuration
    // Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1

    common::Grid3D grid(config.get<IndexType>("NZ"), config.get<IndexType>("NY"), config.get<IndexType>("NX"));
    // distribute the grid onto available processors
    dmemo::DistributionPtr dist(new dmemo::GridDistribution(grid, commShot));

    // Create an object of the mapping (3D-1D) class Coordinates
    Acquisition::Coordinates<ValueType> modelCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    derivatives->init(dist, ctx, config, modelCoordinates, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
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

    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, dist);
    model->prepareForModelling(modelCoordinates, ctx, dist, commShot);
    CheckParameter::checkNumericalArtefeactsAndInstabilities<ValueType>(config, sourceSettings, *model, commAll);

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

    // calculate vector with unique shot numbers and get number of shots
    std::vector<scai::IndexType> uniqueShotNos;
    calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = uniqueShotNos.size();

    /* general block distribution of shot domains accorting to their weights */

    IndexType firstShot = 0;
    IndexType lastShot   = numshots - 1;

    float processorWeight = 1.0f;
    common::Settings::getEnvironment(processorWeight, "WEIGHT");
    float domainWeight = commShot->sum(processorWeight);

    if ( commShot->getRank() == 0)
    {
         // master processors of shot domains determine the load distribution
         auto shotDist = dmemo::genBlockDistributionByWeight(numshots, domainWeight, commInterShot );
         firstShot = shotDist->lb();
         lastShot = shotDist->ub();
    }

    commShot->bcast( &firstShot, 1, 0 );
    commShot->bcast( &lastShot, 1, 0 );

    for (IndexType shotInd = firstShot; shotInd < lastShot; shotInd++) {
        IndexType shotNumber = uniqueShotNos[shotInd];
        /* Update Source */
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
        Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
        sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

        if (config.get<bool>("useReceiversPerShot")) {
            receivers.init(config, modelCoordinates, ctx, dist, shotNumber);
        }

        HOST_PRINT(commShot, "Start time stepping for shot " << shotInd << " (shot no: " << shotNumber << "), domain = " << domain << "\n", 
                             "\nTotal Number of time steps: " << tStepEnd << "\n");
        wavefields->resetWavefields();

        start_t = common::Walltime::get();

        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

            if (tStep % 100 == 0 && tStep != 0) {
                HOST_PRINT(commShot, " ", "Calculating time step " << tStep << " in shot  " << shotNumber << "\n");
            }

            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);

            if (config.get<IndexType>("snapType") > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), DT) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT)) % Common::time2index(config.get<ValueType>("tincSnapshot"), DT) == 0) {
                wavefields->write(config.get<IndexType>("snapType"), config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".", tStep, *derivatives, *model, config.get<IndexType>("PartitionedOut"));
            }
        }

        end_t = common::Walltime::get();
        HOST_PRINT(commShot, "Finished time stepping for shot no: " << shotNumber << " in " << end_t - start_t << " sec.\n");

        receivers.getSeismogramHandler().normalize();

        receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);

        solver->resetCPML();
    }
    return 0;
}
