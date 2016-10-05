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

#include "Configuration.hpp"
#include "Derivatives.hpp"
#include "InitializeMatrices.hpp"

#include "Timesteps.hpp"

#include "Modelparameter/Modelparameter3Dacoustic.hpp"
#include "Wavefields/Wavefields3Dacoustic.hpp"

#include "Receivers.hpp"
#include "Sources.hpp"

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


int main( int argc, char* argv[] )
{
    typedef double ValueType;
    
    if(argc!=2){
        std::cout<< "\n\nNo configuration file given!\n\n" << std::endl;
        return(2);
    }
    
    // read configuration parameter from file
    Configuration<ValueType> config(argv[1]);

    // LAMA specific configuration variables

    // execution context
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    // inter node communicator
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    // inter node distribution
    // block distribution: i-st processor gets lines [i * N/num_processes] to [(i+1) * N/num_processes - 1] of the matrix
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( config.getN(), comm ) );
    
    HOST_PRINT( comm, "\nAcoustic 3-D FD-Algorithm\n\n" );
    if( comm->getRank() == MASTER )
    {
        config.print();
    }

    
    // for timing
    double start_t, end_t;

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    lama::CSRSparseMatrix<ValueType> A, B, C, D, E, F;
    initializeMatrices( A, B, C, D, E, F, dist, ctx, config.getNX(), config.getNY(), config.getNZ(), comm );
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished initializing matrices in " << end_t - start_t << " sec.\n\n" );
    
    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    Wavefields3Dacoustic<ValueType> wavefields(ctx,dist);
    
    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    
    /* Receivers */
    Receivers<ValueType> receivers(config,dist);
    
    /* Sources */
    Sources<ValueType> sources(config,dist);
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    Modelparameter3Dacoustic<ValueType> model(config,ctx,dist);

    /* --------------------------------------- */
    /* Time stepping                           */
    /* --------------------------------------- */
    HOST_PRINT( comm, "Start time stepping\n" );

    start_t = common::Walltime::get();
    
    timesteps( receivers, sources, model, wavefields, A, B, C, D, E, F, config.getVfactor(), config.getPfactor(), config.getNT(), lama::Scalar( 1.0/config.getDH() ), comm, dist );
    
    end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );

    
    receivers.writeSeismograms(config.getSeismogramFilename());
    

    return 0;
}


