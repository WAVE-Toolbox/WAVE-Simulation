#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration/Configuration.hpp"

#include "Modelparameter/Modelparameter3Dacoustic.hpp"
#include "Wavefields/Wavefields3Dacoustic.hpp"

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"

#include "ForwardSolver/Derivatives/Derivatives3D.hpp"

#include "ForwardSolver/ForwardSolver.hpp"

#include "ForwardSolver/ForwardSolver3Dacoustic.hpp"

#include "ForwardSolver/BoundaryConditions/Abs3D.hpp"

#include "Common/HostPrint.hpp"

#include "Partitioning/Partitioning3DCubes.hpp"

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
    Configuration::Configuration<ValueType> config(argv[1]);
    
    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    /* inter node communicator */
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    /* inter node distribution */
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( config.getN(), comm ) );
    
    if( config.getUseCubePartitioning()){
        Partitioning::Partitioning3DCubes<ValueType> partitioning(config,comm);
        dmemo::DistributionPtr dist=partitioning.getDist();
    }
    
    HOST_PRINT( comm, "\nSOFI3D acoustic - LAMA Version\n\n" );
    if( comm->getRank() == MASTER )
    {
        config.print();
    }
    
    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Derivatives::FD3D<ValueType> derivatives( dist, ctx, config, comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n" );
    
    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    ForwardSolver::Abs::Abs3D<ValueType> absorbing_coefficients( dist, ctx, config, comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrix in " << end_t - start_t << " sec.\n\n" );
    
    
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields::FD3Dacoustic<ValueType> wavefields(ctx,dist);
    
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers(config,dist);
    Acquisition::Sources<ValueType> sources(config,dist);
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter::FD3Dacoustic<ValueType> model(config,ctx,dist);
    
    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    
    ForwardSolver::FD3Dacoustic<ValueType> solver;
    
    HOST_PRINT( comm, "Start time stepping\n" );
    start_t = common::Walltime::get();
    
    solver.run( receivers, sources, model, wavefields, derivatives, absorbing_coefficients, config.getNT(), comm);
    
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );
    
    solver.seismogram.writeToFileRaw(config.getSeismogramFilename());
    
    
    return 0;
}


