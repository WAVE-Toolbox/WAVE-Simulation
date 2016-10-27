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
                
                FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, dmemo::CommunicatorPtr comm );
                FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm);
                
                void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, dmemo::CommunicatorPtr comm );
                
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
                
                void derivatives(IndexType NX, IndexType NY, IndexType NZ, IndexType spatialFDorder, dmemo::DistributionPtr dist );
                
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
KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, dmemo::CommunicatorPtr comm ){
    initializeMatrices(dist, ctx, NX, NY, NZ, DH, DT, spatialFDorder, comm );
}


//! \brief Calculation of second order accurate derivatives
/*!
 *
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param dist Distribution of the wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::derivatives(IndexType NX, IndexType NY, IndexType NZ, IndexType spatialFDorder, dmemo::DistributionPtr dist )
{
    SCAI_REGION( "derivatives" )
    
    dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);
    

    /* Do some calculations to avoid it within the for loops */
    IndexType NXNY=NX*NY;
    IndexType N=NX*NY*NZ;
    
    IndexType numLocalIndices=localIndices.size(); //< Number of local indices
    IndexType A_numLocalValues=numLocalIndices; //< Number of local values of Matrix A (here the diagonal elements are added directly)
    IndexType B_numLocalValues=numLocalIndices; //< Number of local values of Matrix B (here the diagonal elements are added directly)
    IndexType C_numLocalValues=numLocalIndices; //< Number of local values of Matrix C (here the diagonal elements are added directly)
    
    /* Add the number of non-diagonal values */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); //< Get read access to localIndices
    IndexType read_localIndices_temp; //< Temporary storage, so we do not have to access the array
    IndexType read_localIndices_temp_plusOne; //< Temporary storage, to save floating point operations
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        read_localIndices_temp_plusOne=read_localIndices_temp+1;
        
        /* Check for non-diagonal elements of A */
        switch (spatialFDorder)
        {
        case 12:
        	if( read_localIndices_temp_plusOne % NX >= 6){	//Check if grid point 5 steps backward is available
        		A_numLocalValues++;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-6) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 6 steps forward is available
        		A_numLocalValues++;
        	}
        case 10:
        	if( read_localIndices_temp_plusOne % NX >= 5){	//Check if grid point 4 steps backward is available
        		A_numLocalValues++;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-5) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 5 steps forward is available
        		A_numLocalValues++;
        	}
        case 8:
        	if( read_localIndices_temp_plusOne % NX >= 4){	//Check if grid point 3 steps backward is available
        		A_numLocalValues++;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-4) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 4 steps forward is available
        		A_numLocalValues++;
        	}
        case 6:
        	if( read_localIndices_temp_plusOne % NX >= 3){	//Check if grid point 2 steps backward is available
        		A_numLocalValues++;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-3) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 3 steps forward is available
        		A_numLocalValues++;
        	}
        case 4:
            if( read_localIndices_temp_plusOne % NX != 1 ){	//Check if grid point 1 step backward is available
                A_numLocalValues++;
            }
            if( (read_localIndices_temp_plusOne % NX <= NX-2) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 2 steps forward is available
            	A_numLocalValues++;
            }
        case 2:
            if( read_localIndices_temp_plusOne % NX != 0 ){	//Check if grid point 1 step forward is available
                A_numLocalValues++;
            }
        }

        
        /* Check for non-diagonal elements of B */
        switch (spatialFDorder)
        {
        case 12:
        	if( ( read_localIndices_temp_plusOne % NXNY > 5*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 5 steps backward is available
        		B_numLocalValues++;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 6) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 6 steps forward is available
        		B_numLocalValues++;
        	}
        case 10:
        	if( ( read_localIndices_temp_plusOne % NXNY > 4*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 4 steps backward is available
        		B_numLocalValues++;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 5) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 5 steps forward is available
        		B_numLocalValues++;
        	}
        case 8:
        	if( ( read_localIndices_temp_plusOne % NXNY > 3*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 3 steps backward is available
        		B_numLocalValues++;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 4) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 4 steps forward is available
        		B_numLocalValues++;
        	}
        case 6:
        	if( ( read_localIndices_temp_plusOne % NXNY > 2*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 2 steps backward is available
        		B_numLocalValues++;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 3) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 3 steps forward is available
        		B_numLocalValues++;
        	}
        case 4:
        	if( ( read_localIndices_temp_plusOne % NXNY > NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 1 step backward is available
        		B_numLocalValues++;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 2) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 2 steps forward is available
        		B_numLocalValues++;
        	}
        case 2:
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 1)) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 1 step forward is available
        		B_numLocalValues++;
        	}
        }
        
        /* Check for non-diagonal elements of C */
        switch (spatialFDorder)
        {
        case 12:
        	if( read_localIndices_temp_plusOne > 5*NXNY ){	//Check if grid point 5 steps backward is available
        		C_numLocalValues++;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 6) ){	//Check if grid point 6 steps forward is available
        		C_numLocalValues++;
        	}
        case 10:
        	if( read_localIndices_temp_plusOne > 4*NXNY ){	//Check if grid point 4 steps backward is available
        		C_numLocalValues++;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 5) ){	//Check if grid point 5 steps forward is available
        		C_numLocalValues++;
        	}
        case 8:
        	if( read_localIndices_temp_plusOne > 3*NXNY ){	//Check if grid point 3 steps backward is available
        		C_numLocalValues++;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 4) ){	//Check if grid point 4 steps forward is available
        		C_numLocalValues++;
        	}
        case 6:
        	if( read_localIndices_temp_plusOne > 2*NXNY ){	//Check if grid point 2 steps backward is available
        		C_numLocalValues++;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 3) ){	//Check if grid point 3 steps forward is available
        		C_numLocalValues++;
        	}
        case 4:
        	if( read_localIndices_temp_plusOne > NXNY ){	//Check if grid point 1 step backward is available
        		C_numLocalValues++;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 2) ){	//Check if grid point 2 steps forward is available
        		C_numLocalValues++;
        	}
        case 2:
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ-1)){	//Check if grid point 1 step forward is available
        		C_numLocalValues++;
        	}
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
    
    ValueType forward_1, forward_2, forward_3, forward_4, forward_5, forward_6;
    ValueType backward_6, backward_5, backward_4, backward_3, backward_2, backward_1;

    // set Taylor FD coefficients for each order
    switch (spatialFDorder)
    {
    case 2:
    	forward_1=1.0;
    	backward_1=-1.0;
    	break;
    case 4:
    	forward_2=-1.0/24.0;
    	forward_1=9.0/8.0;
    	backward_1=-9.0/8.0;
    	backward_2=1.0/24.0;
    	break;
    case 6:
    	forward_3=3.0/640.0;
        forward_2=-25.0/384.0;
       	forward_1=75.0/64.0;
       	backward_1=-75.0/64.0;
       	backward_2=25.0/384.0;
       	backward_3=-3.0/640.0;
       	break;
    case 8:
    	forward_4=-5.0/7168.0;
    	forward_3=49.0/5120.0;
    	forward_2=-245.0/3072.0;
    	forward_1=1225.0/1024.0;
    	backward_1=-1225.0/1024.0;
    	backward_2=245.0/3072.0;
    	backward_3=-49.0/5120.0;
    	backward_4=5.0/7168.0;
    	break;
    case 10:
    	forward_5=8756999275442633.0/73786976294838206464.0;
    	forward_4=-8142668969129685.0/4611686018427387904.0;
    	forward_3=567.0/40960.0;
    	forward_2=-735.0/8192.0;
    	forward_1=19845.0/16384.0;
    	backward_1=-19845.0/16384.0;
    	backward_2=735.0/8192.0;
    	backward_3=-567.0/40960.0;
    	backward_4=8142668969129685.0/4611686018427387904.0;
    	backward_5=-8756999275442633.0/73786976294838206464.0;
    	break;
    case 12:
    	forward_6=-6448335830095439.0/295147905179352825856.0;
    	forward_5=1655620175512543.0/4611686018427387904.0;
    	forward_4=-6842103786556949.0/2305843009213693952.0;
    	forward_3=628618285389933.0/36028797018963968.0;
    	forward_2=-436540475965291.0/4503599627370496.0;
    	forward_1=2750204998582123.0/2251799813685248.0;
    	backward_1=-2750204998582123.0/2251799813685248.0;
    	backward_2=436540475965291.0/4503599627370496.0;
    	backward_3=-628618285389933.0/36028797018963968.0;
    	backward_4=6842103786556949.0/2305843009213693952.0;
    	backward_5=-1655620175512543.0/4611686018427387904.0;
    	backward_6=6448335830095439.0/295147905179352825856.0;
    }

    /* Set the values into the indice arrays and the value array for CSR matrix build */
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        read_localIndices_temp_plusOne=read_localIndices_temp+1;
        
        /*----------*/
        /* Matrix A */
        /*----------*/

        switch (spatialFDorder)
        {
        case 12:
        	if( read_localIndices_temp_plusOne % NX >= 6){	//Check if grid point 5 steps backward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp-5;
        		A_write_valuesLocal[A_countJA]=backward_6;
        		A_countJA++;
        	}
        case 10:
        	if( read_localIndices_temp_plusOne % NX >= 5){	//Check if grid point 4 steps backward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp-4;
        		A_write_valuesLocal[A_countJA]=backward_5;
        		A_countJA++;
       		}
     	case 8:
        	if( read_localIndices_temp_plusOne % NX >= 4){	//Check if grid point 3 steps backward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp-3;
        		A_write_valuesLocal[A_countJA]=backward_4;
        		A_countJA++;
       		}
       	case 6:
           	if( read_localIndices_temp_plusOne % NX >= 3){	//Check if grid point 2 steps backward is available
           		A_write_csrJALocal[A_countJA]=read_localIndices_temp-2;
           		A_write_valuesLocal[A_countJA]=backward_3;
           		A_countJA++;
            }
        case 4:
        	if( read_localIndices_temp_plusOne % NX != 1 ){	//Check if grid point 1 step backward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp-1;
        		A_write_valuesLocal[A_countJA]=backward_2;
        		A_countJA++;
        	}
        case 2:
        	A_write_csrJALocal[A_countJA]=read_localIndices_temp;	//diagonal element
        	A_write_valuesLocal[A_countJA]=backward_1;
        	A_countJA++;
        	if( read_localIndices_temp_plusOne % NX != 0 ){	//Check if grid point 1 step forward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp+1;
        		A_write_valuesLocal[A_countJA]=forward_1;
        		A_countJA++;
        	}
        	if(spatialFDorder==2){
        		break;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-2) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 2 steps forward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp+2;
        		A_write_valuesLocal[A_countJA]=forward_2;
        		A_countJA++;
        	}
        	if(spatialFDorder==4){
        		break;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-3) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 3 steps forward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp+3;
        		A_write_valuesLocal[A_countJA]=forward_3;
        	    A_countJA++;
        	}
        	if(spatialFDorder==6){
        		break;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-4) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 4 steps forward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp+4;
        		A_write_valuesLocal[A_countJA]=forward_4;
        		A_countJA++;
        	}
        	if(spatialFDorder==8){
        		break;
        	}
       		if( (read_localIndices_temp_plusOne % NX <= NX-5) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 5 steps forward is available
       			A_write_csrJALocal[A_countJA]=read_localIndices_temp+5;
       			A_write_valuesLocal[A_countJA]=forward_5;
       			A_countJA++;
       		}
        	if(spatialFDorder==10){
        		break;
        	}
        	if( (read_localIndices_temp_plusOne % NX <= NX-6) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 6 steps forward is available
        		A_write_csrJALocal[A_countJA]=read_localIndices_temp+6;
        		A_write_valuesLocal[A_countJA]=forward_6;
        		A_countJA++;
        	}
        }
        A_write_csrIALocal[A_countIA]=A_countJA;
        A_countIA++;

        /*----------*/
        /* Matrix B */
        /*----------*/
        switch (spatialFDorder)
        {
        case 12:
        	if( ( read_localIndices_temp_plusOne % NXNY > 5*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 5 steps backward is available
                B_write_csrJALocal[B_countJA]=read_localIndices_temp-5*NX;
                B_write_valuesLocal[B_countJA]=backward_6;
                B_countJA++;
        	}
        case 10:
        	if( ( read_localIndices_temp_plusOne % NXNY > 4*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 4 steps backward is available
                B_write_csrJALocal[B_countJA]=read_localIndices_temp-4*NX;
                B_write_valuesLocal[B_countJA]=backward_5;
                B_countJA++;
        	}
        case 8:
        	if( ( read_localIndices_temp_plusOne % NXNY > 3*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 3 steps backward is available
                B_write_csrJALocal[B_countJA]=read_localIndices_temp-3*NX;
                B_write_valuesLocal[B_countJA]=backward_4;
                B_countJA++;
        	}
        case 6:
        	if( ( read_localIndices_temp_plusOne % NXNY > 2*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 2 steps backward is available
                B_write_csrJALocal[B_countJA]=read_localIndices_temp-2*NX;
                B_write_valuesLocal[B_countJA]=backward_3;
                B_countJA++;
        	}
        case 4:
        	if( ( read_localIndices_temp_plusOne %  NXNY  > NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 1 step backward is available
                B_write_csrJALocal[B_countJA]=read_localIndices_temp-NX;
                B_write_valuesLocal[B_countJA]=backward_2;
                B_countJA++;
        	}
        case 2:
        	B_write_csrJALocal[B_countJA]=read_localIndices_temp;	//diagonal element
        	B_write_valuesLocal[B_countJA]=backward_1;
        	B_countJA++;
        	if( ( read_localIndices_temp_plusOne %  NXNY  <= NX * (NY - 1)) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 1 step forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+NX;
        		B_write_valuesLocal[B_countJA]=forward_1;
        		B_countJA++;
        	}
        	if(spatialFDorder==2){
        		break;
        	}
        	if( ( read_localIndices_temp_plusOne %  NXNY  <= NX * (NY - 2) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 2 steps forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+2*NX;
        		B_write_valuesLocal[B_countJA]=forward_2;
        		B_countJA++;
        	}
        	if(spatialFDorder==4){
        		break;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 3) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 3 steps forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+3*NX;
        		B_write_valuesLocal[B_countJA]=forward_3;
        		B_countJA++;
        	}
        	if(spatialFDorder==6){
        		break;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 4) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 4 steps forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+4*NX;
        		B_write_valuesLocal[B_countJA]=forward_4;
        		B_countJA++;
        	}
        	if(spatialFDorder==8){
        		break;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 5) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 5 steps forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+5*NX;
        		B_write_valuesLocal[B_countJA]=forward_5;
        		B_countJA++;
        	}
        	if(spatialFDorder==10){
        		break;
        	}
        	if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 6) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 6 steps forward is available
        		B_write_csrJALocal[B_countJA]=read_localIndices_temp+6*NX;
        		B_write_valuesLocal[B_countJA]=forward_6;
        		B_countJA++;
        	}
        }
        B_write_csrIALocal[B_countIA]=B_countJA;
        B_countIA++;

        /*----------*/
        /* Matrix C */
        /*----------*/

        switch (spatialFDorder)
        {
        case 12:
        	if( read_localIndices_temp_plusOne > 5*NXNY ){	//Check if grid point 5 steps backward is available
                C_write_csrJALocal[C_countJA]=read_localIndices_temp-5*NXNY;
                C_write_valuesLocal[C_countJA]=backward_6;
                C_countJA++;
        	}
        case 10:
        	if( read_localIndices_temp_plusOne > 4*NXNY ){	//Check if grid point 4 steps backward is available
                C_write_csrJALocal[C_countJA]=read_localIndices_temp-4*NXNY;
                C_write_valuesLocal[C_countJA]=backward_5;
                C_countJA++;
        	}
        case 8:
        	if( read_localIndices_temp_plusOne > 3*NXNY ){	//Check if grid point 3 steps backward is available
                C_write_csrJALocal[C_countJA]=read_localIndices_temp-3*NXNY;
                C_write_valuesLocal[C_countJA]=backward_4;
                C_countJA++;
        	}
        case 6:
        	if( read_localIndices_temp_plusOne > 2*NXNY ){	//Check if grid point 2 steps backward is available
                C_write_csrJALocal[C_countJA]=read_localIndices_temp-2*NXNY;
                C_write_valuesLocal[C_countJA]=backward_3;
                C_countJA++;
        	}
        case 4:
        	if( read_localIndices_temp_plusOne > NXNY ){	//Check if grid point 1 step backward is available
                C_write_csrJALocal[C_countJA]=read_localIndices_temp-NXNY;
                C_write_valuesLocal[C_countJA]=backward_2;
                C_countJA++;
        	}
        case 2:
        	C_write_csrJALocal[C_countJA]=read_localIndices_temp;	//diagonal element
        	C_write_valuesLocal[C_countJA]=backward_1;
        	C_countJA++;
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 1) ){	//Check if grid point 1 step forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+NXNY;
        		C_write_valuesLocal[C_countJA]=forward_1;
        		C_countJA++;
        	}
        	if(spatialFDorder==2){
        		break;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 2) ){	//Check if grid point 2 steps forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+2*NXNY;
        		C_write_valuesLocal[C_countJA]=forward_2;
        		C_countJA++;
        	}
        	if(spatialFDorder==4){
        		break;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 3) ){	//Check if grid point 3 steps forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+3*NXNY;
        		C_write_valuesLocal[C_countJA]=forward_3;
        		C_countJA++;
        	}
        	if(spatialFDorder==6){
        		break;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 4) ){	//Check if grid point 4 steps forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+4*NXNY;
        		C_write_valuesLocal[C_countJA]=forward_4;
        		C_countJA++;
        	}
        	if(spatialFDorder==8){
        		break;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 5) ){	//Check if grid point 5 steps forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+5*NXNY;
        		C_write_valuesLocal[C_countJA]=forward_5;
        		C_countJA++;
        	}
        	if(spatialFDorder==10){
        		break;
        	}
        	if( read_localIndices_temp_plusOne <= NXNY * (NZ - 6) ){	//Check if grid point 6 steps forward is available
        		C_write_csrJALocal[C_countJA]=read_localIndices_temp+6*NXNY;
        		C_write_valuesLocal[C_countJA]=forward_6;
        		C_countJA++;
        	}
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
    lama::CSRStorage<ValueType> A_LocalCSR(numLocalIndices,N,A_numLocalValues,A_csrIALocal,A_csrJALocal,A_valuesLocal);
    A.assign(A_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix B, than create distributed CSR matrix B */
    lama::CSRStorage<ValueType> B_LocalCSR(numLocalIndices,N,B_numLocalValues,B_csrIALocal,B_csrJALocal,B_valuesLocal);
    B.assign(B_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix C, than create distributed CSR matrix C */
    lama::CSRStorage<ValueType> C_LocalCSR(numLocalIndices,N,C_numLocalValues,C_csrIALocal,C_csrJALocal,C_valuesLocal);
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
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, dmemo::CommunicatorPtr comm )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the matrices A,B,C,D,E,F...\n" );
    
    derivatives(NX, NY, NZ, spatialFDorder, dist);
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
