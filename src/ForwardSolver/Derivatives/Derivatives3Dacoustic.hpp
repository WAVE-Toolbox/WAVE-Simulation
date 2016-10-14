#pragma once

#include "Derivatives.hpp"
#include "../../Common/HostPrint.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
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
                
                lama::CSRSparseMatrix<ValueType>& getA();
                lama::CSRSparseMatrix<ValueType>& getB();
                lama::CSRSparseMatrix<ValueType>& getC();
                lama::CSRSparseMatrix<ValueType>& getD();
                lama::CSRSparseMatrix<ValueType>& getE();
                lama::CSRSparseMatrix<ValueType>& getF();
                
            private:
                
                lama::CSRSparseMatrix<ValueType> A; //!< Derivative matrix A
                lama::CSRSparseMatrix<ValueType> B; //!< Derivative matrix B
                lama::CSRSparseMatrix<ValueType> C; //!< Derivative matrix C
                lama::CSRSparseMatrix<ValueType> D; //!< Derivative matrix D
                lama::CSRSparseMatrix<ValueType> E; //!< Derivative matrix E
                lama::CSRSparseMatrix<ValueType> F; //!< Derivative matrix F
                
                void derivatives(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist );
                
            };
        } /* end namespace Derivatives */
    } /* end namespace ForwardSolver */
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
KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
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
KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm ){
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
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::derivatives(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist )
{
    SCAI_REGION( "derivatives" )
    
    
    dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);
    
    IndexType numLocalIndices=localIndices.size(); //< Number of local indices
    IndexType A_numLocalValues=numLocalIndices; //< Number of local values of Matrix A (here the diagonal elements are added directly)
    IndexType B_numLocalValues=numLocalIndices; //< Number of local values of Matrix B (here the diagonal elements are added directly)
    IndexType C_numLocalValues=numLocalIndices; //< Number of local values of Matrix C (here the diagonal elements are added directly)
    
    /* Add the number of non-diagonal values */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); //< Get read access to localIndices
    IndexType read_localIndices_temp; //< Temporary storage, so we do not have to access the array
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        
        /* Check for non-diagonal elements of A */
        if( (read_localIndices_temp+1) % NX != 0 ){
            A_numLocalValues++;
        }
        
        /* Check for non-diagonal elements of B */
        if( (( (read_localIndices_temp+1) %  (NX*NY) ) + NX - 1 < NX*NY) &&  ((read_localIndices_temp+1) % (NX*NY) != 0) ){
            B_numLocalValues++;
        }
        
        /* Check for non-diagonal elements of C */
        if( ( read_localIndices_temp+1)  + NX*NY - 1 < NX*NY*NZ &&  (read_localIndices_temp+1) % (NX*NY*NZ) != 0 ){
            C_numLocalValues++;
        }
        
    }
    
    /* Allocate local part to create local CSR storage*/
    hmemo::HArray<ValueType> A_valuesLocal(A_numLocalValues);
    hmemo::HArray<IndexType> A_csrJALocal(A_numLocalValues);
    hmemo::HArray<IndexType> A_csrIALocal(numLocalIndices+1);
    
    hmemo::HArray<ValueType> B_valuesLocal(B_numLocalValues);
    hmemo::HArray<IndexType> B_csrJALocal(B_numLocalValues);
    hmemo::HArray<IndexType> B_csrIALocal(numLocalIndices+1);
    
    hmemo::HArray<ValueType> C_valuesLocal(C_numLocalValues);
    hmemo::HArray<IndexType> C_csrJALocal(C_numLocalValues);
    hmemo::HArray<IndexType> C_csrIALocal(numLocalIndices+1);
    
    /* Get WriteAccess to local part */
    hmemo::WriteAccess<IndexType> A_write_csrJALocal(A_csrJALocal);
    hmemo::WriteAccess<IndexType> A_write_csrIALocal(A_csrIALocal);
    hmemo::WriteAccess<ValueType> A_write_valuesLocal(A_valuesLocal);
    
    hmemo::WriteAccess<IndexType> B_write_csrJALocal(B_csrJALocal);
    hmemo::WriteAccess<IndexType> B_write_csrIALocal(B_csrIALocal);
    hmemo::WriteAccess<ValueType> B_write_valuesLocal(B_valuesLocal);
    
    hmemo::WriteAccess<IndexType> C_write_csrJALocal(C_csrJALocal);
    hmemo::WriteAccess<IndexType> C_write_csrIALocal(C_csrIALocal);
    hmemo::WriteAccess<ValueType> C_write_valuesLocal(C_valuesLocal);
    
    /* Set some counters to create the CSR Storage */
    IndexType A_countJA=0;
    IndexType A_countIA=0;
    A_write_csrIALocal[0]=0;
    A_countIA++;
    
    IndexType B_countJA=0;
    IndexType B_countIA=0;
    B_write_csrIALocal[0]=0;
    B_countIA++;
    
    IndexType C_countJA=0;
    IndexType C_countIA=0;
    C_write_csrIALocal[0]=0;
    C_countIA++;
    
    /* Set the values into the indice arrays and the value array */
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        
        /*----------*/
        /* Matrix A */
        /*----------*/
        
        /* Add diagonal element */
        A_write_csrJALocal[A_countJA]=read_localIndices_temp;
        A_write_valuesLocal[A_countJA]=-1;
        A_countJA++;
        
        /* Add non-diagonal element */
        if( (read_localIndices_temp+1) % NX != 0 ){
            A_write_csrJALocal[A_countJA]=read_localIndices_temp+1;
            A_write_valuesLocal[A_countJA]=1;
            A_countJA++;
        }
        A_write_csrIALocal[A_countIA]=A_countJA;
        A_countIA++;
        
        
        /*----------*/
        /* Matrix B */
        /*----------*/
        
        /* Add diagonal element */
        B_write_csrJALocal[B_countJA]=read_localIndices_temp;
        B_write_valuesLocal[B_countJA]=-1;
        B_countJA++;
        
        /* Add non-diagonal element */
        if( ( (read_localIndices_temp+1) %  (NX*NY) ) + NX - 1 < NX*NY &&  (read_localIndices_temp+1) % (NX*NY) != 0 ){            B_write_csrJALocal[B_countJA]=read_localIndices_temp+NX;
            B_write_valuesLocal[B_countJA]=1;
            B_countJA++;
        }
        B_write_csrIALocal[B_countIA]=B_countJA;
        B_countIA++;
        
        /*----------*/
        /* Matrix C */
        /*----------*/
        
        /* Add diagonal element */
        C_write_csrJALocal[C_countJA]=read_localIndices_temp;
        C_write_valuesLocal[C_countJA]=-1;
        C_countJA++;
        
        /* Add non-diagonal element */
        if( ( read_localIndices_temp+1)  + NX*NY - 1 < NX*NY*NZ &&  (read_localIndices_temp+1) % (NX*NY*NZ) != 0 ){
            C_write_csrJALocal[C_countJA]=read_localIndices_temp+NX*NY;
            C_write_valuesLocal[C_countJA]=1;
            C_countJA++;
        }
        C_write_csrIALocal[C_countIA]=C_countJA;
        C_countIA++;
        
        
    }
    
    /* Release all read and write access */
    read_localIndices.release();
    
    A_write_csrJALocal.release();
    A_write_csrIALocal.release();
    A_write_valuesLocal.release();
    
    B_write_csrJALocal.release();
    B_write_csrIALocal.release();
    B_write_valuesLocal.release();
    
    C_write_csrJALocal.release();
    C_write_csrIALocal.release();
    C_write_valuesLocal.release();
    
    /* Create local CSR storage of Matrix A, than create distributed CSR matrix A */
    lama::CSRStorage<ValueType> A_LocalCSR(numLocalIndices,NX*NY*NZ,A_numLocalValues,A_csrIALocal,A_csrJALocal,A_valuesLocal);
    A.assign(A_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix B, than create distributed CSR matrix B */
    lama::CSRStorage<ValueType> B_LocalCSR(numLocalIndices,NX*NY*NZ,B_numLocalValues,B_csrIALocal,B_csrJALocal,B_valuesLocal);
    B.assign(B_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix C, than create distributed CSR matrix C */
    lama::CSRStorage<ValueType> C_LocalCSR(numLocalIndices,NX*NY*NZ,C_numLocalValues,C_csrIALocal,C_csrJALocal,C_valuesLocal);
    C.assign(C_LocalCSR,dist,dist);
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
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the matrices A,B,C,D,E,F...\n" );
    
    derivatives(NX, NY, NZ, dist);
    HOST_PRINT( comm, "Matrix A, B and C finished.\n" );

    A.setContextPtr( ctx );
    B.setContextPtr( ctx );
    C.setContextPtr( ctx );
    D.setContextPtr( ctx );
    E.setContextPtr( ctx );
    F.setContextPtr( ctx );
    
    D.assignTranspose( A );
    D.scale( -1.0 );
    
    E.assignTranspose( B );
    E.scale( -1.0 );
    
    F.assignTranspose( C );
    F.scale( -1.0 );
    
    HOST_PRINT( comm, "Matrix D, E and F finished.\n" );
    
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
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getA(){
    return(A);
}

//! \brief Getter method for derivative matrix B
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getB(){
    return(B);
}

//! \brief Getter method for derivative matrix C
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getC(){
    return(C);
}

//! \brief Getter method for derivative matrix D
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getD(){
    return(D);
}

//! \brief Getter method for derivative matrix E
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getE(){
    return(E);
}

//! \brief Getter method for derivative matrix F
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::getF(){
    return(F);
}
