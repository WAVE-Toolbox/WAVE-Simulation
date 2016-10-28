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
                
                FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm );
                FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm);
                
                void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm );
                
            private:
                
                using Derivatives<ValueType>::Dxf;	//Dxf: derivative in the direction of the x axis (f: forward) see Documentation for detailed explanation
                using Derivatives<ValueType>::Dyf;
                using Derivatives<ValueType>::Dzf;
                using Derivatives<ValueType>::Dxb;	//b: backward
                using Derivatives<ValueType>::Dyb;
                using Derivatives<ValueType>::Dzb;
                
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
 \param DT Temporal sampling interval#
 \param spatialFDorderInput FD-order of spatial derivative stencils
 \param comm Communicator
 */
template<typename ValueType>
KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::FD3D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm ){
    initializeMatrices(dist, ctx, NX, NY, NZ, DH, DT, spatialFDorderInput, comm );
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
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::derivatives(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist )
{
    SCAI_REGION( "derivatives" )
    
    dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);		//here the local indices of each process are retrieved and stored in localIndices
    
    
    /* Do some calculations to avoid it within the for loops */
    IndexType NXNY=NX*NY;
    IndexType N=NX*NY*NZ;
    
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    IndexType Dxf_numLocalValues=numLocalIndices; // Number of local values of Matrix Dxf (here the diagonal elements are added directly)
    IndexType Dyf_numLocalValues=numLocalIndices; // Number of local values of Matrix Dyf (here the diagonal elements are added directly)
    IndexType Dzf_numLocalValues=numLocalIndices; // Number of local values of Matrix Dzf (here the diagonal elements are added directly)
    
    /* Add the number of non-diagonal values */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp; // Temporary storage of the local index for the ongoing iterations
    IndexType read_localIndices_temp_plusOne; // shifted, so that the definition of the local index starts at 1, instead of 0
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        read_localIndices_temp_plusOne=read_localIndices_temp+1;
        
        /* Check for non-diagonal elements of Dxf */
        switch (this->spatialFDorder)				//here the fall through characteristic of the switch case statement is exploited to run all conditions underneath a certain case
        {
            case 12:
                if( read_localIndices_temp_plusOne % NX >= 6){	//Check if grid point 5 steps backward is available
                    Dxf_numLocalValues++;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-6) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 6 steps forward is available
                    Dxf_numLocalValues++;
                }
            case 10:
                if( read_localIndices_temp_plusOne % NX >= 5){	//Check if grid point 4 steps backward is available
                    Dxf_numLocalValues++;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-5) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 5 steps forward is available
                    Dxf_numLocalValues++;
                }
            case 8:
                if( read_localIndices_temp_plusOne % NX >= 4){	//Check if grid point 3 steps backward is available
                    Dxf_numLocalValues++;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-4) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 4 steps forward is available
                    Dxf_numLocalValues++;
                }
            case 6:
                if( read_localIndices_temp_plusOne % NX >= 3){	//Check if grid point 2 steps backward is available
                    Dxf_numLocalValues++;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-3) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 3 steps forward is available
                    Dxf_numLocalValues++;
                }
            case 4:
                if( read_localIndices_temp_plusOne % NX != 1 ){	//Check if grid point 1 step backward is available
                    Dxf_numLocalValues++;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-2) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 2 steps forward is available
                    Dxf_numLocalValues++;
                }
            case 2:
                if( read_localIndices_temp_plusOne % NX != 0 ){	//Check if grid point 1 step forward is available
                    Dxf_numLocalValues++;
                }
        }
        
        
        /* Check for non-diagonal elements of Dyf */
        switch (this->spatialFDorder)
        {
            case 12:
                if( ( read_localIndices_temp_plusOne % NXNY > 5*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 5 steps backward is available
                    Dyf_numLocalValues++;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 6) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 6 steps forward is available
                    Dyf_numLocalValues++;
                }
            case 10:
                if( ( read_localIndices_temp_plusOne % NXNY > 4*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 4 steps backward is available
                    Dyf_numLocalValues++;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 5) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 5 steps forward is available
                    Dyf_numLocalValues++;
                }
            case 8:
                if( ( read_localIndices_temp_plusOne % NXNY > 3*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 3 steps backward is available
                    Dyf_numLocalValues++;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 4) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 4 steps forward is available
                    Dyf_numLocalValues++;
                }
            case 6:
                if( ( read_localIndices_temp_plusOne % NXNY > 2*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 2 steps backward is available
                    Dyf_numLocalValues++;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 3) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 3 steps forward is available
                    Dyf_numLocalValues++;
                }
            case 4:
                if( ( read_localIndices_temp_plusOne % NXNY > NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 1 step backward is available
                    Dyf_numLocalValues++;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 2) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 2 steps forward is available
                    Dyf_numLocalValues++;
                }
            case 2:
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 1)) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 1 step forward is available
                    Dyf_numLocalValues++;
                }
        }
        
        /* Check for non-diagonal elements of Dzf */
        switch (this->spatialFDorder)
        {
            case 12:
                if( read_localIndices_temp_plusOne > 5*NXNY ){	//Check if grid point 5 steps backward is available
                    Dzf_numLocalValues++;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 6) ){	//Check if grid point 6 steps forward is available
                    Dzf_numLocalValues++;
                }
            case 10:
                if( read_localIndices_temp_plusOne > 4*NXNY ){	//Check if grid point 4 steps backward is available
                    Dzf_numLocalValues++;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 5) ){	//Check if grid point 5 steps forward is available
                    Dzf_numLocalValues++;
                }
            case 8:
                if( read_localIndices_temp_plusOne > 3*NXNY ){	//Check if grid point 3 steps backward is available
                    Dzf_numLocalValues++;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 4) ){	//Check if grid point 4 steps forward is available
                    Dzf_numLocalValues++;
                }
            case 6:
                if( read_localIndices_temp_plusOne > 2*NXNY ){	//Check if grid point 2 steps backward is available
                    Dzf_numLocalValues++;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 3) ){	//Check if grid point 3 steps forward is available
                    Dzf_numLocalValues++;
                }
            case 4:
                if( read_localIndices_temp_plusOne > NXNY ){	//Check if grid point 1 step backward is available
                    Dzf_numLocalValues++;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 2) ){	//Check if grid point 2 steps forward is available
                    Dzf_numLocalValues++;
                }
            case 2:
                if( read_localIndices_temp_plusOne <= NXNY * (NZ-1)){	//Check if grid point 1 step forward is available
                    Dzf_numLocalValues++;
                }
        }
    }
    /* Allocate local part to create local CSR storage*/
    hmemo::HArray<ValueType> Dxf_valuesLocal(Dxf_numLocalValues);
    hmemo::HArray<IndexType> Dxf_csrJALocal(Dxf_numLocalValues);
    hmemo::HArray<IndexType> Dxf_csrIALocal(numLocalIndices+1);
    
    hmemo::HArray<ValueType> Dyf_valuesLocal(Dyf_numLocalValues);
    hmemo::HArray<IndexType> Dyf_csrJALocal(Dyf_numLocalValues);
    hmemo::HArray<IndexType> Dyf_csrIALocal(numLocalIndices+1);
    
    hmemo::HArray<ValueType> Dzf_valuesLocal(Dzf_numLocalValues);
    hmemo::HArray<IndexType> Dzf_csrJALocal(Dzf_numLocalValues);
    hmemo::HArray<IndexType> Dzf_csrIALocal(numLocalIndices+1);
    
    /* Get WriteAccess to local part */
    hmemo::WriteAccess<IndexType> Dxf_write_csrJALocal(Dxf_csrJALocal);
    hmemo::WriteAccess<IndexType> Dxf_write_csrIALocal(Dxf_csrIALocal);
    hmemo::WriteAccess<ValueType> Dxf_write_valuesLocal(Dxf_valuesLocal);
    
    hmemo::WriteAccess<IndexType> Dyf_write_csrJALocal(Dyf_csrJALocal);
    hmemo::WriteAccess<IndexType> Dyf_write_csrIALocal(Dyf_csrIALocal);
    hmemo::WriteAccess<ValueType> Dyf_write_valuesLocal(Dyf_valuesLocal);
    
    hmemo::WriteAccess<IndexType> Dzf_write_csrJALocal(Dzf_csrJALocal);
    hmemo::WriteAccess<IndexType> Dzf_write_csrIALocal(Dzf_csrIALocal);
    hmemo::WriteAccess<ValueType> Dzf_write_valuesLocal(Dzf_valuesLocal);
    
    /* Set some counters to create the CSR Storage */
    IndexType Dxf_countJA=0;
    IndexType Dxf_countIA=0;
    Dxf_write_csrIALocal[0]=0;
    Dxf_countIA++;
    
    IndexType Dyf_countJA=0;
    IndexType Dyf_countIA=0;
    Dyf_write_csrIALocal[0]=0;
    Dyf_countIA++;
    
    IndexType Dzf_countJA=0;
    IndexType Dzf_countIA=0;
    Dzf_write_csrIALocal[0]=0;
    Dzf_countIA++;
    
    ValueType forward_1, forward_2, forward_3, forward_4, forward_5, forward_6;
    ValueType backward_6, backward_5, backward_4, backward_3, backward_2, backward_1;
    
    // set Taylor FD coefficients for each order
    switch (this->spatialFDorder)
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
            break;
        default:
            COMMON_THROWEXCEPTION(" Unkown spatialFDorder value.")
            break;
    }
    
    /* Set the values into the indice arrays and the value array for CSR matrix build */
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        read_localIndices_temp_plusOne=read_localIndices_temp+1;
        
        /*------------*/
        /* Matrix Dxf */
        /*------------*/
        
        switch (this->spatialFDorder)
        {
            case 12:
                if( read_localIndices_temp_plusOne % NX >= 6){	//Check if grid point 5 steps backward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp-5;
                    Dxf_write_valuesLocal[Dxf_countJA]=backward_6;
                    Dxf_countJA++;
                }
            case 10:
                if( read_localIndices_temp_plusOne % NX >= 5){	//Check if grid point 4 steps backward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp-4;
                    Dxf_write_valuesLocal[Dxf_countJA]=backward_5;
                    Dxf_countJA++;
                }
            case 8:
                if( read_localIndices_temp_plusOne % NX >= 4){	//Check if grid point 3 steps backward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp-3;
                    Dxf_write_valuesLocal[Dxf_countJA]=backward_4;
                    Dxf_countJA++;
                }
            case 6:
                if( read_localIndices_temp_plusOne % NX >= 3){	//Check if grid point 2 steps backward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp-2;
                    Dxf_write_valuesLocal[Dxf_countJA]=backward_3;
                    Dxf_countJA++;
                }
            case 4:
                if( read_localIndices_temp_plusOne % NX != 1 ){	//Check if grid point 1 step backward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp-1;
                    Dxf_write_valuesLocal[Dxf_countJA]=backward_2;
                    Dxf_countJA++;
                }
            case 2:
                Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp;	//diagonal element
                Dxf_write_valuesLocal[Dxf_countJA]=backward_1;
                Dxf_countJA++;
                if( read_localIndices_temp_plusOne % NX != 0 ){	//Check if grid point 1 step forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+1;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_1;
                    Dxf_countJA++;
                }
                if(this->spatialFDorder==2){
                    break;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-2) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 2 steps forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+2;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_2;
                    Dxf_countJA++;
                }
                if(this->spatialFDorder==4){
                    break;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-3) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 3 steps forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+3;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_3;
                    Dxf_countJA++;
                }
                if(this->spatialFDorder==6){
                    break;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-4) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 4 steps forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+4;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_4;
                    Dxf_countJA++;
                }
                if(this->spatialFDorder==8){
                    break;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-5) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 5 steps forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+5;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_5;
                    Dxf_countJA++;
                }
                if(this->spatialFDorder==10){
                    break;
                }
                if( (read_localIndices_temp_plusOne % NX <= NX-6) && (read_localIndices_temp_plusOne % NX != 0)){	//Check if grid point 6 steps forward is available
                    Dxf_write_csrJALocal[Dxf_countJA]=read_localIndices_temp+6;
                    Dxf_write_valuesLocal[Dxf_countJA]=forward_6;
                    Dxf_countJA++;
                }
        }
        Dxf_write_csrIALocal[Dxf_countIA]=Dxf_countJA;
        Dxf_countIA++;
        
        /*------------*/
        /* Matrix Dyf */
        /*------------*/
        switch (this->spatialFDorder)
        {
            case 12:
                if( ( read_localIndices_temp_plusOne % NXNY > 5*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 5 steps backward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp-5*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=backward_6;
                    Dyf_countJA++;
                }
            case 10:
                if( ( read_localIndices_temp_plusOne % NXNY > 4*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 4 steps backward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp-4*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=backward_5;
                    Dyf_countJA++;
                }
            case 8:
                if( ( read_localIndices_temp_plusOne % NXNY > 3*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 3 steps backward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp-3*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=backward_4;
                    Dyf_countJA++;
                }
            case 6:
                if( ( read_localIndices_temp_plusOne % NXNY > 2*NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 2 steps backward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp-2*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=backward_3;
                    Dyf_countJA++;
                }
            case 4:
                if( ( read_localIndices_temp_plusOne %  NXNY  > NX ) || (read_localIndices_temp_plusOne % NXNY == 0) ){	//Check if grid point 1 step backward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp-NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=backward_2;
                    Dyf_countJA++;
                }
            case 2:
                Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp;	//diagonal element
                Dyf_write_valuesLocal[Dyf_countJA]=backward_1;
                Dyf_countJA++;
                if( ( read_localIndices_temp_plusOne %  NXNY  <= NX * (NY - 1)) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 1 step forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_1;
                    Dyf_countJA++;
                }
                if(this->spatialFDorder==2){
                    break;
                }
                if( ( read_localIndices_temp_plusOne %  NXNY  <= NX * (NY - 2) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 2 steps forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+2*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_2;
                    Dyf_countJA++;
                }
                if(this->spatialFDorder==4){
                    break;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 3) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 3 steps forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+3*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_3;
                    Dyf_countJA++;
                }
                if(this->spatialFDorder==6){
                    break;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 4) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 4 steps forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+4*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_4;
                    Dyf_countJA++;
                }
                if(this->spatialFDorder==8){
                    break;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 5) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 5 steps forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+5*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_5;
                    Dyf_countJA++;
                }
                if(this->spatialFDorder==10){
                    break;
                }
                if( ( read_localIndices_temp_plusOne % NXNY <= NX * (NY - 6) ) && (read_localIndices_temp_plusOne % NXNY != 0) ){	//Check if grid point 6 steps forward is available
                    Dyf_write_csrJALocal[Dyf_countJA]=read_localIndices_temp+6*NX;
                    Dyf_write_valuesLocal[Dyf_countJA]=forward_6;
                    Dyf_countJA++;
                }
        }
        Dyf_write_csrIALocal[Dyf_countIA]=Dyf_countJA;
        Dyf_countIA++;
        
        /*------------*/
        /* Matrix Dzf */
        /*------------*/
        
        switch (this->spatialFDorder)
        {
            case 12:
                if( read_localIndices_temp_plusOne > 5*NXNY ){	//Check if grid point 5 steps backward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp-5*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=backward_6;
                    Dzf_countJA++;
                }
            case 10:
                if( read_localIndices_temp_plusOne > 4*NXNY ){	//Check if grid point 4 steps backward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp-4*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=backward_5;
                    Dzf_countJA++;
                }
            case 8:
                if( read_localIndices_temp_plusOne > 3*NXNY ){	//Check if grid point 3 steps backward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp-3*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=backward_4;
                    Dzf_countJA++;
                }
            case 6:
                if( read_localIndices_temp_plusOne > 2*NXNY ){	//Check if grid point 2 steps backward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp-2*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=backward_3;
                    Dzf_countJA++;
                }
            case 4:
                if( read_localIndices_temp_plusOne > NXNY ){	//Check if grid point 1 step backward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp-NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=backward_2;
                    Dzf_countJA++;
                }
            case 2:
                Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp;	//diagonal element
                Dzf_write_valuesLocal[Dzf_countJA]=backward_1;
                Dzf_countJA++;
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 1) ){	//Check if grid point 1 step forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_1;
                    Dzf_countJA++;
                }
                if(this->spatialFDorder==2){
                    break;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 2) ){	//Check if grid point 2 steps forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+2*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_2;
                    Dzf_countJA++;
                }
                if(this->spatialFDorder==4){
                    break;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 3) ){	//Check if grid point 3 steps forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+3*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_3;
                    Dzf_countJA++;
                }
                if(this->spatialFDorder==6){
                    break;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 4) ){	//Check if grid point 4 steps forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+4*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_4;
                    Dzf_countJA++;
                }
                if(this->spatialFDorder==8){
                    break;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 5) ){	//Check if grid point 5 steps forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+5*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_5;
                    Dzf_countJA++;
                }
                if(this->spatialFDorder==10){
                    break;
                }
                if( read_localIndices_temp_plusOne <= NXNY * (NZ - 6) ){	//Check if grid point 6 steps forward is available
                    Dzf_write_csrJALocal[Dzf_countJA]=read_localIndices_temp+6*NXNY;
                    Dzf_write_valuesLocal[Dzf_countJA]=forward_6;
                    Dzf_countJA++;
                }
        }
        Dzf_write_csrIALocal[Dzf_countIA]=Dzf_countJA;
        Dzf_countIA++;
    }
    
    /* Release all read and write access */
    read_localIndices.release();
    
    Dxf_write_csrJALocal.release();
    Dxf_write_csrIALocal.release();
    Dxf_write_valuesLocal.release();
    
    Dyf_write_csrJALocal.release();
    Dyf_write_csrIALocal.release();
    Dyf_write_valuesLocal.release();
    
    Dzf_write_csrJALocal.release();
    Dzf_write_csrIALocal.release();
    Dzf_write_valuesLocal.release();
    
    /* Create local CSR storage of Matrix Dxf, than create distributed CSR matrix Dxf */
    lama::CSRStorage<ValueType> Dxf_LocalCSR(numLocalIndices,N,Dxf_numLocalValues,Dxf_csrIALocal,Dxf_csrJALocal,Dxf_valuesLocal);
    Dxf.assign(Dxf_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix Dyf, than create distributed CSR matrix Dyf */
    lama::CSRStorage<ValueType> Dyf_LocalCSR(numLocalIndices,N,Dyf_numLocalValues,Dyf_csrIALocal,Dyf_csrJALocal,Dyf_valuesLocal);
    Dyf.assign(Dyf_LocalCSR,dist,dist);
    
    /* Create local CSR storage of Matrix Dzf, than create distributed CSR matrix Dzf */
    lama::CSRStorage<ValueType> Dzf_LocalCSR(numLocalIndices,N,Dzf_numLocalValues,Dzf_csrIALocal,Dzf_csrJALocal,Dzf_valuesLocal);
    Dzf.assign(Dzf_LocalCSR,dist,dist);
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
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FD3D<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the matrices Dxf, Dyf, Dzf, Dxb, Dyb, Dzb…\n" );
    
    // Set FD-order to class member
    this->spatialFDorder=spatialFDorderInput;
    
    derivatives(NX, NY, NZ, dist);
    HOST_PRINT( comm, "Matrix Dxf, Dyf and Dzf finished.\n" );
    
    Dxf.setContextPtr( ctx );
    Dyf.setContextPtr( ctx );
    Dzf.setContextPtr( ctx );
    Dxb.setContextPtr( ctx );
    Dyb.setContextPtr( ctx );
    Dzb.setContextPtr( ctx );
    
    Dxb.assignTranspose( Dxf );
    Dxb.scale( -1.0 );
    
    Dyb.assignTranspose( Dyf );
    Dyb.scale( -1.0 );
    
    Dzb.assignTranspose( Dzf );
    Dzb.scale( -1.0 );
    
    HOST_PRINT( comm, "Matrix Dxb, Dyb and Dzb finished.\n" );
    
    Dxf.scale(lama::Scalar(DT/DH));
    Dyf.scale(lama::Scalar(DT/DH));
    Dzf.scale(lama::Scalar(DT/DH));
    Dxb.scale(lama::Scalar(DT/DH));
    Dyb.scale(lama::Scalar(DT/DH));
    Dzb.scale(lama::Scalar(DT/DH));
    
    HOST_PRINT( comm, "Finished with initialization of the matrices!\n" );
    
    
}

