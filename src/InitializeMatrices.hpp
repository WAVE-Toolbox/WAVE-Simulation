
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

/*
 *  routine for initializing all system matrices
 *  uses derivatives for initializing A, B, C
 */
template<typename ValueType>
void initializeMatrices( lama::SparseMatrix<ValueType>& A, lama::SparseMatrix<ValueType>& B, lama::SparseMatrix<ValueType>& C,
                        lama::Matrix& D, lama::Matrix& E, lama::Matrix& F, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,
                        IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the matrices A,B,C,D,E,F...\n" );
    
    derivatives( A, B, C, NX, NY, NZ, dist, ctx, comm );
    
    D.setContextPtr( ctx );
    E.setContextPtr( ctx );
    F.setContextPtr( ctx );
    
    D.assignTranspose( A );
    D.scale( -1.0 );
    HOST_PRINT( comm, "Matrix D finished\n" );
    
    E.assignTranspose( B );
    E.scale( -1.0 );
    HOST_PRINT( comm, "Matrix E finished\n" );
    
    F.assignTranspose( C );
    F.scale( -1.0 );
    HOST_PRINT( comm, "Matrix F finished\n" );
    
    A.scale(lama::Scalar(DT/DH));
    B.scale(lama::Scalar(DT/DH));
    C.scale(lama::Scalar(DT/DH));
    D.scale(lama::Scalar(DT/DH));
    E.scale(lama::Scalar(DT/DH));
    F.scale(lama::Scalar(DT/DH));
    
    HOST_PRINT( comm, "Finished with initialization of the matrices!\n" );
}


