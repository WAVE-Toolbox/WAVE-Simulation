#pragma once

#include "Derivatives.hpp"

#define MASTER 0

#define HOST_PRINT( comm, msg )     \
{                                   \
int myRank = comm->getRank();   \
if ( myRank == MASTER )         \
{                               \
std::cout << msg;           \
}                               \
}

namespace KITGPI {
    
    //! \brief Derivatives namespace
    namespace Derivatives {
        
        //! \brief Class for the calculation of the Derivatives matrices for 3-D FD Simulations
        /*!
         * Calculation of derivative matrices for an equidistand grid
         *
         */
        template<typename ValueType>
        class FD3D : public Derivatives<ValueType>
        {
        public:
            
            //! Default constructor
            FD3D(){};
            
            //! Default destructor
            ~FD3D(){};
            
            FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm );
            FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm);
            
            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm );

            lama::CSRSparseMatrix<ValueType>* getA();
            lama::CSRSparseMatrix<ValueType>* getB(); 
            lama::CSRSparseMatrix<ValueType>* getC(); 
            lama::CSRSparseMatrix<ValueType>* getD();
            lama::CSRSparseMatrix<ValueType>* getE();
            lama::CSRSparseMatrix<ValueType>* getF();
            
        private:
            
            lama::CSRSparseMatrix<ValueType> A; //!< Derivative matrix A
            lama::CSRSparseMatrix<ValueType> B; //!< Derivative matrix B
            lama::CSRSparseMatrix<ValueType> C; //!< Derivative matrix C
            lama::CSRSparseMatrix<ValueType> D; //!< Derivative matrix D
            lama::CSRSparseMatrix<ValueType> E; //!< Derivative matrix E
            lama::CSRSparseMatrix<ValueType> F; //!< Derivative matrix F
            
            void derivatives(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, dmemo::CommunicatorPtr comm );
            
        };
    } /* end namespace Derivatives */
} /* end namespace KITGPI */


//! \brief Constructor to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
KITGPI::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    Derivatives<ValueType>::initializeMatrices(dist,ctx, config, comm );
}

//! \brief Constructor of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval
 \param comm Communicator
 */
template<typename ValueType>
KITGPI::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm ){
    initializeMatrices(dist, ctx, NX, NY,NZ,DH, DT, comm );
}


//! \brief Calculation of second order accurate derivatives
/*!
 *
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::Derivatives::FD3D<ValueType>::derivatives(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, dmemo::CommunicatorPtr comm )
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
    
    /* --------------- */
    /* create matrix A */
    /* --------------- */
    numValues = NX + (NX - 1); // diagonal element (NX) + secondary diagonal elements (NX - 1)
    csrIA.reserve( NX + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );
    
    IndexType count = 0;
    IndexType size = NX;
    csrIA.push_back( 0 );
    for( IndexType i = 0; i < size; ++i ){
        if(i<size){
            csrJA.push_back(i);
            csrValues.push_back( -1.0 );
            count++;
        }
        if((i+1)<size){
            csrJA.push_back(i+1);
            csrValues.push_back( +1.0 );
            count++;
        }
        csrIA.push_back( count );
    }
    
    /* create CSR storage help */
    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    
    /* deallocate memory of std::vectors, in order to free memory for the next operation  */
    csrIA = std::vector<IndexType>();
    csrJA = std::vector<IndexType>();
    csrValues = std::vector<ValueType>();
    
    /* create matrix A */
    lama::MatrixCreator::buildReplicatedDiag( A, *storageHelp, NZ * NY ); // not distributed, need to redistribute afterwards
    A.redistribute( dist, dist );
    A.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix A finished\n" );
    
    /* deallocate storageHelp, in order to free memory for the next operation */
    storageHelp->purge();

    /* --------------- */
    /* create matrix B */
    /* --------------- */
    numValues = NX*NY + (NX*NY - (NX+1) + 1); // diagonal element (NZ*NX) + secondary diagonal elements (NZ*NX - (NZ+1) + 1)
    csrIA.reserve( NX + 1 );
    csrJA.reserve( numValues );
    csrValues.reserve( numValues );
    
    count = 0;
    size = NX * NY;
    csrIA.push_back( 0 );
    for( IndexType i = 0; i < size; ++i ){
        if(i<size){
            csrJA.push_back(i);
            csrValues.push_back( -1.0 );
            count++;
        }
        if(i+NX<size){
            csrJA.push_back(i+NX);
            csrValues.push_back( +1.0 );
            count++;
        }
        csrIA.push_back( count );
    }
    
    storageHelp->setRawCSRData( size, size, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );
    
    /* deallocate memory of std::vectors, in order to free memory for the next operation  */
    csrIA = std::vector<IndexType>();
    csrJA = std::vector<IndexType>();
    csrValues = std::vector<ValueType>();
    
    lama::MatrixCreator::buildReplicatedDiag( B, *storageHelp, NZ ); // not distributed, need to redistribute afterwards
    B.redistribute( dist, dist );
    B.setContextPtr( ctx );
    HOST_PRINT( comm, "Matrix B finished\n" );
    
    /* deallocate storageHelp, in order to free memory for the next operation */
    storageHelp->purge();
    
    /* --------------- */
    /* create matrix C */
    /* --------------- */
    // initialize by diagonals
    int myRank   = comm->getRank();
    int numRanks = comm->getSize();
    
    IndexType globalSize = dist->getGlobalSize();
    IndexType numDiagonals = 2;
    IndexType secondaryIndex = NX * NY; // = remaining part of secondary diagonal
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

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::Derivatives::FD3D<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the matrices A,B,C,D,E,F...\n" );
    
    derivatives(NX, NY, NZ, dist, ctx, comm );
    
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

//! \brief Getter method for derivative matrix A
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getA(){
    return(&A);
}

//! \brief Getter method for derivative matrix B
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getB(){
    return(&B);
}

//! \brief Getter method for derivative matrix C
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getC(){
    return(&C);
}

//! \brief Getter method for derivative matrix D
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getD(){
    return(&D);
}

//! \brief Getter method for derivative matrix E
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getE(){
    return(&E);
}

//! \brief Getter method for derivative matrix F
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>* KITGPI::Derivatives::FD3D<ValueType>::getF(){
    return(&F);
}
