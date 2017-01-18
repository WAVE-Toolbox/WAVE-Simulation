#include <scai/lama.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration/Configuration.hpp"

#include "Modelparameter/Elastic.hpp"
#include "Wavefields/Wavefields2Delastic.hpp"

#include "Acquisition/Sources.hpp"
#include "Acquisition/Receivers.hpp"

#include "ForwardSolver/ForwardSolver.hpp"
#include "ForwardSolver/ForwardSolver2Delastic.hpp"

#include "ForwardSolver/Derivatives/FDTD2D.hpp"
#include "ForwardSolver/BoundaryCondition/FreeSurface3Delastic.hpp"

#include "Common/HostPrint.hpp"
#include "Partitioning/PartitioningCubes.hpp"

using namespace scai;
using namespace KITGPI;

int main( int argc, char* argv[] )
{
    typedef double ValueType;
    double start_t, end_t; /* For timing */
    
    if(argc!=2){
        std::cout<< "\n\nNo configuration file given!\n\n" << std::endl;
        return(2);
    }
    
    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);
    
    SCAI_ASSERT( config.get<IndexType>("NZ") == 1 , " NZ must be equal to 1 for 2-D simulation" );
    
    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */ 
    /* inter node communicator */
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank( comm->getNodeRank() );
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    /* inter node distribution */
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    IndexType getN = config.get<IndexType>("NZ") * config.get<IndexType>("NX") * config.get<IndexType>("NY");
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( getN, comm ) );
    
    if( config.get<IndexType>("UseCubePartitioning")){
        Partitioning::PartitioningCubes<ValueType> partitioning(config,comm);
        dist=partitioning.getDist();
    }
    
    HOST_PRINT( comm, "\nSOFI2D elastic - LAMA Version\n\n" );
    if( comm->getRank() == MASTER )
    {
        config.print();
    }
    
    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::FDTD2D<ValueType> derivatives( dist, ctx, config, comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n" );
    
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::FD2Delastic<ValueType> wavefields(ctx,dist);
    
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers(config,ctx,dist);
    Acquisition::Sources<ValueType> sources(config,ctx,dist);
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::Elastic<ValueType> model(config,ctx,dist);
    
    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    
    ForwardSolver::FD2Delastic<ValueType> solver;
    
    solver.prepareBoundaryConditions(config,derivatives,dist,ctx);
    
    IndexType getNT = static_cast<IndexType>( ( config.get<ValueType>("T") / config.get<ValueType>("DT") ) + 0.5 );
    solver.run(receivers, sources, model, wavefields, derivatives, getNT, config.get<ValueType>("DT"));
    
    receivers.getSeismogramHandler().write(config);

    return 0;
}


