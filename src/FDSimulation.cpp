#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath> 

#include "Configuration.hpp"
#include "Derivatives.hpp"
#include "InitializeMatrices.hpp"
#include "SourceFunction.hpp"
#include "Timesteps.hpp"
#include "Modelparameter.hpp"

#include "Wavefields/Wavefields3Dacoustic.hpp"

using namespace scai;

#define MASTER 0

#define HOST_PRINT( comm, msg )     \
{                                   \
int myRank = comm->getRank();   \
if ( myRank == MASTER )         \
{                               \
std::cout << msg;           \
}                               \
}


int main( int /*argc*/, char** /*argv[]*/ )
{
    // we do all calculation in double precision
    typedef double ValueType;

    // read configuration parameter from file
    Configuration<ValueType> config( "input/Configuration.txt" );

    // LAMA specific configuration variables

    // execution context
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    // inter node communicator
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    // inter node distribution
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( config.getN(), comm ) );

    HOST_PRINT( comm, "Acoustic 3-D FD-Algorithm\n\n" );
    if( comm->getRank() == MASTER )
    {
        config.print();
    }

    
    // for timing
    double start_t, end_t;
   
    /* --------------------------------------- */
    /* Calculate source signal                 */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    lama::DenseVector<ValueType> source( config.getNT(), ValueType(0), config.getDT(), ctx );
    sourceFunction( source, config.getFC(), config.getAMP(), comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished calculating source in " << end_t - start_t << " sec.\n\n" );

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    lama::CSRSparseMatrix<ValueType> A, B, C, D, E, F;
    initializeMatrices( A, B, C, D, E, F, dist, ctx, config.getNX(), config.getNY(), config.getNZ(), comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n" );
    
    
    Wavefields3Dacoustic<ValueType> wavefields(ctx,dist);
    
    // seismogram data: to store at each time step
    lama::DenseVector<ValueType> seismogram( config.getNT(), 0.0 ); // no ctx, use default: Host
    // Model
    Modelparameter<ValueType> model;
    if(config.getReadModel()==1){
        model.init(ctx,dist,config.getFilenameModel());
    } else {
        model.init(ctx,dist,config.getM(),config.getRho());
    }
    
    HOST_PRINT( comm, "Start time stepping\n" );

    /* --------------------------------------- */
    /* Time stepping                           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    timesteps( seismogram, source, model, wavefields, A, B, C, D, E, F,
               config.getVfactor(), config.getPfactor(), config.getNT(), lama::Scalar( 1.0/config.getDH() ),
               config.getSourceIndex(), config.getSeismogramIndex(), comm, dist );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );

    // print vector data for seismogram plot
    seismogram.writeToFile( "seismograms/seismogram.mtx" );


    return 0;
}


