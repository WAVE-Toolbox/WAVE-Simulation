#include <scai/lama.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/GridDistribution.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration/Configuration.hpp"

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"

#include "ForwardSolver/ForwardSolver.hpp"

#include "ForwardSolver/Derivatives/DerivativesFactory.hpp"
#include "ForwardSolver/ForwardSolverFactory.hpp"
#include "Modelparameter/ModelparameterFactory.hpp"
#include "Wavefields/WavefieldsFactory.hpp"

#include "Common/HostPrint.hpp"
#include "Partitioning/PartitioningCubes.hpp"

#include "CheckParameter/CheckParameter.hpp"

using namespace scai;
using namespace KITGPI;

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA
 
    common::Settings::parseArgs( argc, argv );
 
    typedef float ValueType;
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
    
    IndexType firstSnapshot = 0, lastSnapshot = 0, incSnapshot = 0 ;
    if (config.get<IndexType>("snapType") > 0){
        firstSnapshot = static_cast<IndexType>(config.get<ValueType>("tFirstSnapshot") / config.get<ValueType>("DT") + 0.5);
        lastSnapshot = static_cast<IndexType>(config.get<ValueType>("tLastSnapshot") / config.get<ValueType>("DT") + 0.5);
        incSnapshot = static_cast<IndexType>(config.get<ValueType>("tIncSnapshot") / config.get<ValueType>("DT") + 0.5);
    }
    
    IndexType partitionedOut = static_cast<IndexType>(config.get<ValueType>("PartitionedOut"));

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* inter node communicator */
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(comm->getNodeRank());
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    // following lines should be part of checkParameter.tpp
    if (comm->getSize() != config.get<IndexType>("ProcNX") * config.get<IndexType>("ProcNY") * config.get<IndexType>("ProcNZ")) {
        HOST_PRINT(comm, "\n Error: Number of MPI processes (" << comm->getSize() << ") doesn't match the number of processes specified in " << argv[1] << ": ProcNX*ProcNY*ProcNZ=" << config.get<IndexType>("ProcNX")*config.get<IndexType>("ProcNY")*config.get<IndexType>("ProcNZ") << "\n")
        return(2);	    
    }
    
    // inter node distribution 
    // define the grid topology by sizes NX, NY, and NZ from configuration
    // Attention: LAMA uses row-major indexing while SOFI-3D uses column-major, so switch dimensions, x-dimension has stride 1

    common::Grid3D grid(config.get<IndexType>("NZ"), config.get<IndexType>("NY"), config.get<IndexType>("NX") );
    common::Grid3D procGrid(config.get<IndexType>("ProcNZ"), config.get<IndexType>("ProcNY"), config.get<IndexType>("ProcNX") );
    // distribute the grid onto available processors, topology can be set by environment variable
    dmemo::DistributionPtr dist(new dmemo::GridDistribution( grid, comm, procGrid));

    HOST_PRINT(comm, "\nSOFI" << dimension << " " << equationType << " - LAMA Version\n\n");
    if (comm->getRank() == MASTERGPI) {
        config.print();
    }
    
    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    derivatives->init(dist, ctx, config, comm);
    end_t = common::Walltime::get();
    HOST_PRINT(comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    wavefields->init(ctx, dist);

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Sources<ValueType> sources(config, ctx, dist);
    Acquisition::Receivers<ValueType> receivers;
    if (!config.get<bool>("useReceiversPerShot"))
        receivers.init(config, ctx, dist);
    CheckParameter::checkAcquisitionGeometry<ValueType,IndexType>(config,comm,sources.getNumShots());  
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    model->init(config, ctx, dist);
    model->prepareForModelling(config, ctx, dist, comm);
    
    CheckParameter::checkNumericalArtefeactsAndInstabilities<ValueType,IndexType>(config, *model,comm);
    
    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    
    
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    solver->initForwardSolver(config, *derivatives, *wavefields, *model, ctx, config.get<ValueType>("DT"));
    solver->prepareForModelling(*model, config.get<ValueType>("DT"));

    for (IndexType shotNumber = 0; shotNumber < sources.getNumShots(); shotNumber++) {
        /* Update Source */
        if (!config.get<bool>("runSimultaneousShots"))
            sources.init(config, ctx, dist, shotNumber);
        if (config.get<bool>("useReceiversPerShot"))
            receivers.init(config, ctx, dist, shotNumber);

        HOST_PRINT(comm, "Start time stepping for shot " << shotNumber + 1 << " of " << sources.getNumShots() << "\n"
                                                         << "Total Number of time steps: " << tStepEnd << "\n");
        wavefields->resetWavefields();
	
        start_t = common::Walltime::get();
	
        for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {
            
            if (tStep % 100 == 0 && tStep != 0) {
                HOST_PRINT(comm, "Calculating time step " << tStep << "\n");
            }
            
            solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
            
            if (config.get<IndexType>("snapType") > 0 && tStep >= firstSnapshot && tStep <= lastSnapshot && (tStep-firstSnapshot)%incSnapshot == 0) {
                    wavefields->write(config.get<IndexType>("snapType"),config.get<std::string>("WavefieldFileName"),tStep, *derivatives, *model, partitionedOut);
            }
            
        }
        
        end_t = common::Walltime::get();
        HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

        receivers.getSeismogramHandler().normalize();
        
        if (!config.get<bool>("runSimultaneousShots")) {
        receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber));
        } else {
            receivers.getSeismogramHandler().write(config, config.get<std::string>("SeismogramFilename"));
        }
        
        solver->resetCPML();
    }
    return 0;
}

