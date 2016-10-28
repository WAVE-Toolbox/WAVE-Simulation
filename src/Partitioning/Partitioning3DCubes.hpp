#pragma once

#include "Partitioning.hpp"
#include "../Acquisition/Coordinates.hpp"
#include "../Configuration/Configuration.hpp"

using namespace scai;

namespace KITGPI {
    
    namespace Partitioning {
        
        //! \brief Creating a partition of 3-D cubes
        /*!
         * This class can create a partition of the 3-D wavefield that consists of 3-D cubes.
         *
         * procNX*procNY*procNZ have to be equal to comm->getSize();
         *
         */
        template <typename ValueType>
        class Partitioning3DCubes : public Partitioning<ValueType>, protected Acquisition::Coordinates<ValueType>
        {
            
        public:
            
            //! Default constructor
            Partitioning3DCubes(){};
            
            Partitioning3DCubes(Configuration::Configuration<ValueType> config,dmemo::CommunicatorPtr comm);
            
            //! Default destructor
            ~Partitioning3DCubes(){};
            
            dmemo::DistributionPtr getDist();
            
        private:
            
            dmemo::DistributionPtr calculate(IndexType procNX,IndexType procNY,IndexType procNZ,IndexType NX,IndexType NY,IndexType NZ, dmemo::CommunicatorPtr comm);
            
            dmemo::DistributionPtr dist_cubes; //!< Distribution
            
        };
        
    }
}


/*! \brief Getter method for the disttribution pointer *
 *
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Partitioning::Partitioning3DCubes<ValueType>::getDist(){
    if(dist_cubes==NULL){
        COMMON_THROWEXCEPTION ( "Distribution ist not set " )
    }
    return(dist_cubes);
}

/*! \brief Constructor based on the configuration class and the communicator
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::Partitioning::Partitioning3DCubes<ValueType>::Partitioning3DCubes(Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm){
    dist_cubes=calculate(config.getProcNX(),config.getProcNY(),config.getProcNZ(),config.getNX(),config.getNY(),config.getNZ(),comm);
}


//! \brief Calculation of 3-D cube partition
/*!
 * This class can create a partition of the 3-D wavefield that consists of 3-D cubes.
 *
 * procNX*procNY*procNZ have to be equal to comm->getSize();
 *
 \param procNX Number of cores in X-direction
 \param procNY Number of cores in Y-direction
 \param procNZ Number of cores in Z-direction
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in X-direction
 \param NZ Number of grid points in X-direction
 \param comm Communicator for the distribution
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Partitioning::Partitioning3DCubes<ValueType>::calculate(IndexType procNX,IndexType procNY,IndexType procNZ,IndexType NX,IndexType NY,IndexType NZ, dmemo::CommunicatorPtr comm){
    
    
    IndexType rank=comm->getRank();
    IndexType numRanks=comm->getSize();
    
    /* Check some things */
    if( procNX*procNY*procNZ != numRanks ){ COMMON_THROWEXCEPTION(" Number of cores differ from config to actual setting  ") }
    if(NX % procNX != 0 ){ COMMON_THROWEXCEPTION(" NX % procNX != 0  ") }
    if(NY % procNY != 0 ){ COMMON_THROWEXCEPTION(" NY % procNY != 0  ") }
    if(NZ % procNZ != 0 ){ COMMON_THROWEXCEPTION(" NZ % procNZ != 0  ") }
    
    /* Calculate 3-D coordinates of the CPU */
    IndexType posNX= rank % procNX;
    IndexType posNZ= rank / ( procNX * procNY);
    IndexType posNY= (( rank - ( procNX * procNY ) * posNZ ) / procNX );
    
    /* Calculate range of grid points */
    IndexType range_x_lower=( NX / procNX ) * posNX;
    IndexType range_y_lower=( NY / procNY ) * posNY;
    IndexType range_z_lower=( NZ / procNZ ) * posNZ;
    
    IndexType range_x_upper=( NX / procNX ) * (posNX+1);
    IndexType range_y_upper=( NY / procNY ) * (posNY+1);
    IndexType range_z_upper=( NZ / procNZ ) * (posNZ+1);
    
    IndexType numGlobalGridPoints= NX * NY * NZ ;
    IndexType numLocalGridPoints=( numGlobalGridPoints ) / numRanks;
    
    /* Determine local indices */
    hmemo::HArray<IndexType> localIndices;
    localIndices.resize(numLocalGridPoints);
    hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
    IndexType i=0;
    IndexType indice;
    for(IndexType x=0; x<NX; x++){
        for(IndexType y=0; y<NY; y++){
            for(IndexType z=0; z<NZ; z++){
                if( x>=range_x_lower && x<range_x_upper && y>=range_y_lower && y<range_y_upper && z>=range_z_lower && z<range_z_upper ){
                    indice=this->coordinate2index(x+1,y+1,z+1,NX,NY, NZ);
                    write_localIndices[i]=indice;
                    i++;
                }
            }
        }
    }
    if( i != numLocalGridPoints ) {  COMMON_THROWEXCEPTION(" i != numLocalGridPoints   " << i << numLocalGridPoints ); }
    write_localIndices.release();
    
    /* create GeneralDistribution */
    dmemo::DistributionPtr dist_cubus( new dmemo::GeneralDistribution(numGlobalGridPoints,localIndices,comm));
    
    return(dist_cubus);
}
