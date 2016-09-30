
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
 *  routine for initializing derivative matrices (needed for velocity Updates vz, vx, vy)
 *  create csr data structure (sub)matrices (and replicate along matrix diagonal)
 */
template<typename ValueType>
void derivatives( lama::SparseMatrix<ValueType>& A,
                 lama::SparseMatrix<ValueType>& B,
                 lama::SparseMatrix<ValueType>& C,
                 IndexType NX,
                 IndexType NY,
                 IndexType NZ,
                 dmemo::DistributionPtr dist,
                 hmemo::ContextPtr ctx,
                 dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "derivatives" )
    
    // Matrix A,B are created in 2 steps:
    //   1: create MatrixStorage in CSR format
    //   2: duplicate MatrixStorage along Diagonal till full Matrix size is reached
    
    // for creating CSR MatrixStorage
    common::unique_ptr<lama::MatrixStorage<ValueType> > storageHelp( new lama::CSRStorage<ValueType>() );
    IndexType numValues;
    std::vector<IndexType> csrIA;
    std::vector<IndexType> csrJA;
    std::vector<ValueType> csrValues;
    
    // create matrix A
    numValues = NZ + (NZ - 1); // diagonal element (NZ) + secondary diagonal elements (NZ - 1)
    csrIA.reserve( NZ + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );
    
    IndexType count = 0;
    IndexType size = NZ;
    csrIA.push_back( 0 );
    for( IndexType i = 0; i < size; ++i )
    {
        for( IndexType j = 0; j < size; ++j )
        {
            if ( i == j )
            {
                ++count;
                csrJA.push_back(j);
                csrValues.push_back( -1.0 );
            }
            
            if ( j - 1 == i )
            {
                ++count;
                csrJA.push_back(j);
                csrValues.push_back( 1.0 );
            }
        }
        csrIA.push_back( count );
    }
    
    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    lama::MatrixCreator::buildReplicatedDiag( A, *storageHelp, NX * NY ); // not distributed, need to redistribute afterwards
    A.redistribute( dist, dist );
    A.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix A finished\n" );
    
    csrIA.clear();
    csrJA.clear();
    csrValues.clear();
    
    // create matrix B
    numValues = NZ*NX + (NZ*NX - (NZ+1) + 1); // diagonal element (NZ*NX) + secondary diagonal elements (NZ*NX - (NZ+1) + 1)
    csrIA.reserve( NZ + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );
    
    count = 0;
    size = NZ * NX;
    csrIA.push_back( 0 );
    for( IndexType i = 0; i < size; ++i )
    {
        for( IndexType j = 0; j < size; ++j )
        {
            if ( i == j )
            {
                ++count;
                csrJA.push_back(j);
                csrValues.push_back( -1.0 );
            }
            
            if ( j - NZ == i )
            {
                ++count;
                csrJA.push_back(j);
                csrValues.push_back( 1.0 );
            }
        }
        csrIA.push_back( count );
    }
    
    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    lama::MatrixCreator::buildReplicatedDiag( B, *storageHelp, NY ); // not distributed, need to redistribute afterwards
    B.redistribute( dist, dist );
    B.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix B finished\n" );
    
    csrIA.clear();
    csrJA.clear();
    csrValues.clear();
    
    // create matrix C
    // initialize by diagonals
    
    int myRank   = comm->getRank();
    int numRanks = comm->getSize();
    
    IndexType globalSize = dist->getGlobalSize();
    IndexType numDiagonals = 2;
    IndexType secondaryIndex = NZ * NX; // = remaining part of secondary diagonal
    IndexType numSecondary = globalSize - secondaryIndex;
    
    int lb, ub;
    dmemo::BlockDistribution::getLocalRange( lb, ub, dist->getGlobalSize(), myRank, numRanks );
    
    size = dist->getLocalSize(); //getGlobalSize();
    
    IndexType myStart = lb;
    IndexType myEnd   = std::min( ub, numSecondary );
    IndexType mySize  = std::max( myEnd - myStart, 0 );
    IndexType myRemaining = size - mySize;
    
    std::vector<IndexType> offsets;
    // distributed offset start with their lb
    offsets.push_back( lb + 0 );              // index main diagonal (starting with '0')
    offsets.push_back( lb + secondaryIndex ); // index secondary diagonal
    
    std::vector<ValueType> diagonals;
    diagonals.reserve( numDiagonals * size );
    // insert backwards so we can insert from the beginning of the vector
    diagonals.insert( diagonals.begin(), myRemaining, 0.0 ); // add ghost secondary diagonal (out of column bound)
    diagonals.insert( diagonals.begin(), mySize,      1.0 ); // add secondary diagonal
    diagonals.insert( diagonals.begin(), size,       -1.0 ); // insert main diagonal before secondary diagonal
    
    C.setRawDIAData( dist, dist, numDiagonals, &offsets[0], &diagonals[0] );
    C.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix C finished\n" );
}

