#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Coordinates.hpp"


namespace KITGPI {
    
    namespace Acquisition {
        
        //! Handling of receivers
        /*!
         * This class accounts for the handling of seismic receivers.
         * It provides the reading of the receivers acquisition from file, the distribution of the receivers and the collection of the seismograms.
         */
        template <typename ValueType>
        class Receivers : private Coordinates<ValueType>
        {
            
        public:
            
            Receivers(){};
            Receivers(Configuration::Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield);
            ~Receivers(){};
            
            void readReceiverAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield);
            void saveSeismogramsLocal(lama::DenseVector<ValueType>& wavefield, IndexType nt, IndexType NT);
            
            void writeSeismograms(std::string filename);
            void writeReceiverAcquisition(std::string filename);
            
            IndexType getNumReceiversGlobal();
            IndexType getNumReceiversLocal();
            
        private:
            
            void getLocalReceivers(dmemo::DistributionPtr dist_wavefield);
            void getReceiverDistribution(dmemo::CommunicatorPtr comm);
            void allocateSeismograms(IndexType NT);
            
            IndexType numReceiversGlobal=0; //!< Number of receivers global
            IndexType numReceiversLocal=0; //!< Number of receivers local
            
            lama::DenseVector<ValueType> coordinates; //!< Coordinates of receivers global (1-D coordinates)
            
            dmemo::DistributionPtr dist_wavefield_receivers=NULL; //!< Calculated Distribution of the receivers based on the distribution of the wavefields
            dmemo::DistributionPtr no_dist_NT=NULL; //!< No distribution of the columns of the seismogram matrix
            
            //! Seismograms
            lama::DenseMatrix<ValueType> seismograms;
            
            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition; //!< Matrix that stores the receiver acquisition
            IndexType numParameter=0; //!< Number of receiver parameters given in acquisition matrix
            lama::DenseVector<ValueType> receiver_type; //!< Type of Receivers: 1==Pressure
            
        };
    }
}

/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
KITGPI::Acquisition::Receivers<ValueType>::Receivers(Configuration::Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield){
    readReceiverAcquisition(config.getReceiverFilename(),config.getNX(), config.getNY(), config.getNZ(),dist_wavefield);
    allocateSeismograms(config.getNT());
}


/*! \brief Get number of global receivers
 *
 \return Number of global receivers
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversGlobal(){
    return(numReceiversGlobal);
}


/*! \brief Get number of local receivers
 *
 \return Number of local receivers on this process
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversLocal(){
    return(numReceiversLocal);
}


/*! \brief Write receivers acquisition to file
 *
 \param filename Filename to write receivers acquisition
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::writeReceiverAcquisition(std::string filename){
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.assignTranspose(acquisition);
    acquisition_temp.writeToFile(filename);
}


/*! \brief Read receivers acquisition from file
 *
 * This method reads the receivers acquisition, calculates the 1-D coordinates from the 3-D coordinates,
 * splits up the receivers configuration into the corresponding vectors and calculates the receivers distribution.
 * The parameter vectors will be distributed accordingly to the receiver distribution.
 *
 \param filename Filename to read receivers acquisition
 \param NX Number of global grid points in X
 \param NY Number of global grid points in Y
 \param NZ Number of global grid points in Z
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::readReceiverAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield)
{
    
    /* Read acquisition matrix */
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(filename);
    
    IndexType nrow_temp=acquisition_temp.getNumRows();
    IndexType ncolumn_temp=acquisition_temp.getNumColumns();
    
    /* Derive number of receivers and number of read-in parameters */
    numReceiversGlobal=nrow_temp;
    numParameter=ncolumn_temp;
    
    /* Check if number of parameters is supported */
    if( numParameter!=4 ) {
        COMMON_THROWEXCEPTION ( "Receivers acquisition file has an unkown format " )
    }
    
    /* Distribution: Master process only (comm->myRank()==0) */
    dmemo::DistributionPtr dist_master_numParameter( new dmemo::CyclicDistribution( numParameter, numParameter, dist_wavefield->getCommunicatorPtr() ) );
    dmemo::DistributionPtr dist_master_numReceiversGlobal( new dmemo::CyclicDistribution( numReceiversGlobal, numReceiversGlobal, dist_wavefield->getCommunicatorPtr() )  );
    
    /* Distribution: Replicated on all processes */
    dmemo::DistributionPtr no_dist_numReceiversGlobal( new scai::dmemo::NoDistribution ( numReceiversGlobal ) );
    dmemo::DistributionPtr no_dist_numParameter( new scai::dmemo::NoDistribution ( numParameter ) );
    
    /* Allocate acquisition matrix on master */
    acquisition.allocate(dist_master_numParameter,no_dist_numReceiversGlobal);
    
    /* Allocate coordinates on master */
    coordinates.allocate(dist_master_numReceiversGlobal);
    
    /* Allocate receiver parameter vectors on master */
    receiver_type.allocate(dist_master_numReceiversGlobal);
    
    /* Local operations on master: 1. Transpose acquisition, 2. calculate 1-D coordinates  */
    if(dist_wavefield->getCommunicator().getRank()==0){
        
        /* Get WriteAccess to local data of acquisition */
        lama::DenseStorage<ValueType>* acquisition_DS=&acquisition.getLocalStorage();
        hmemo::HArray<ValueType>* acquisition_HA=&acquisition_DS->getData();
        hmemo::WriteAccess<ValueType> write_acquisition_HA(*acquisition_HA);
        
        /* Get Readaccess to local data of acquisition_temp */
        lama::DenseStorage<ValueType>* acquisition_temp_DS=&acquisition_temp.getLocalStorage();
        hmemo::HArray<ValueType>* acquisition_temp_HA=&acquisition_temp_DS->getData();
        hmemo::ReadAccess<ValueType> read_acquisition_temp_HA(*acquisition_temp_HA);
        
        /* Transpose local data */
        for(IndexType row=0; row<nrow_temp;row++){
            for(IndexType column=0; column<ncolumn_temp;column++){
                write_acquisition_HA[ row + nrow_temp*column ] = read_acquisition_temp_HA[ column + ncolumn_temp*row ];
            }
        }
        
        /* Release write and read access to local data */
        write_acquisition_HA.release();
        read_acquisition_temp_HA.release();
        
        /* Get readAccess to acquisition matrix (local) */
        acquisition_DS=&acquisition.getLocalStorage();
        acquisition_HA=&acquisition_DS->getData();
        hmemo::ReadAccess<ValueType> read_acquisition_HA(*acquisition_HA);
        
        /* Get writeAccess to coordinates vector (local) */
        utilskernel::LArray<ValueType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<ValueType> write_coordinates_LA(*coordinates_LA);
        
        /* 2. Calculate 1-D coordinates form 3-D coordinates */
        IndexType X,Y,Z;
        for(IndexType i=0; i<numReceiversGlobal; i++){
            
            X=read_acquisition_HA[ i + numReceiversGlobal*0 ];
            Y=read_acquisition_HA[ i + numReceiversGlobal*1 ];
            Z=read_acquisition_HA[ i + numReceiversGlobal*2 ];
            
            write_coordinates_LA[i]=this->coordinate2index(X,Y,Z,NX,NY,NZ);
        }
        
        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
        
    }
    
    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numReceiversGlobal);
    
    /* Get local receivers from global receivers */
    getLocalReceivers(dist_wavefield);
    
    /* Get receiver distribution */
    getReceiverDistribution(dist_wavefield->getCommunicatorPtr());
    
    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_numParameter,no_dist_numReceiversGlobal);
    
    /* Allocate receiver parameter vectors on all processes */
    receiver_type.allocate(numReceiversGlobal);
    
    /* Save receiver configurations from acquisition matrix in vectors */
    acquisition.getLocalRow(receiver_type,3);
    
    /* Redistribute receiver parameter vectors to corresponding processes */
    coordinates.redistribute(dist_wavefield_receivers);
    receiver_type.redistribute(dist_wavefield_receivers);
    
    
}


/*! \brief Write seismograms signals to file
 *
 \param filename Filename to write seismograms
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::writeSeismograms(std::string filename){
    seismograms.writeToFile(filename);
}


/*! \brief Allocation of the seismogram matrix
 *
 * Allocation of the seismogram matrix based on an already defined receiver distribution and the number of time steps.
 * The seismogram matrix is allocated based on the distributions.
 *
 \param NT Number of time steps
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::allocateSeismograms(IndexType NT)
{
    if(dist_wavefield_receivers==NULL) {
        COMMON_THROWEXCEPTION ( "Row distribution of seismograms (dist_wavefield_receivers) is not set!" )
    }
    
    dmemo::DistributionPtr no_dist_NT_temp( new scai::dmemo::NoDistribution ( NT ) );
    no_dist_NT=no_dist_NT_temp;
    
    /* seismogram matrix is row distributed according to dist_wavefield_receivers, No column distribution */
    seismograms.allocate(dist_wavefield_receivers,no_dist_NT);
}


/*! \brief Calculation of the receiver distribution
 *
 * Calculation of the receiver distribution based on the global number of receivers numReceiversGlobal and the local number of receivers numReceiversLocal on each node.
 * Generates a GenBlockDistribution.
 *
 \param comm Communicator for the generated distribution
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getReceiverDistribution(dmemo::CommunicatorPtr comm)
{
    if(numReceiversGlobal==0){
        COMMON_THROWEXCEPTION ( " There is no global receiver (numReceiversGlobal==0)! ")
    }
    
    dmemo::DistributionPtr dist_temp( new dmemo::GenBlockDistribution(numReceiversGlobal,numReceiversLocal,comm));
    dist_wavefield_receivers=dist_temp;
}


/*! \brief Determine the number of local receivers from the global receivers coordinates
 *
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getLocalReceivers(dmemo::DistributionPtr dist_wavefield)
{
    
    if(coordinates.size()==0){
        COMMON_THROWEXCEPTION ( " The vector coordinates does not contain any elements ! ")
    }
    
    lama::DenseVector<ValueType> coordinateslocal; //!< Coordinates of receivers local (1-D coordinates)
    
    this->Global2Local(coordinates,coordinateslocal,dist_wavefield);
    
    numReceiversLocal=coordinateslocal.size();
    numReceiversGlobal=coordinates.size();
    
}


/*! \brief Records the seismograms during forward modelling (locally)
 *
 * Saves the samples at the specific time.
 *
 \param wavefield Wavefield
 \param nt Current time step
 \param NT Total number of time steps
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::saveSeismogramsLocal(lama::DenseVector<ValueType>& wavefield, IndexType nt, IndexType NT)
{
    if(numReceiversLocal>0){
        
        utilskernel::LArray<ValueType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<ValueType> read_coordinates_LA(*coordinates_LA);
        
        lama::DenseStorage<ValueType>* seismograms_DS=&seismograms.getLocalStorage();
        hmemo::HArray<ValueType>* seismograms_HA=&seismograms_DS->getData();
        hmemo::WriteAccess<ValueType> write_seismograms_HA(*seismograms_HA);
        
        dmemo::DistributionPtr dist_wavefield=wavefield.getDistributionPtr();
        
        scai::lama::Scalar receiver_index_temp;
        IndexType coordinate_global;
        IndexType coordinate_local;
        
        for(IndexType i=0; i<numReceiversLocal; i++){
            coordinate_global=read_coordinates_LA[i];
            coordinate_local=dist_wavefield->global2local(coordinate_global);
            
            write_seismograms_HA[nt+NT*i]=wavefield.getLocalValues()[coordinate_local];
        }
        
        read_coordinates_LA.release();
        write_seismograms_HA.release();
    }
}



