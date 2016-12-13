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
        class Receivers : protected Coordinates<ValueType>
        {
            
        public:
            
            Receivers():numReceiversGlobal(0),numReceiversLocal(0),numParameter(0){};
            explicit Receivers(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);
            ~Receivers(){};
            
            void init(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);
            
            void readReceiverAcquisition(std::string const& filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx);
            void writeReceiverAcquisition(std::string const& filename) const;
            
            /* Getter functions */
            IndexType getNumReceiversGlobal() const;
            IndexType getNumReceiversLocal() const;
            dmemo::DistributionPtr getReceiversDistribution() const;
            lama::DenseVector<IndexType> const& getCoordinates() const;
            lama::DenseVector<IndexType> const& getReceiversType() const;
            
        private:
            
            void getLocalReceivers(dmemo::DistributionPtr dist_wavefield);
            void getReceiverDistribution(dmemo::CommunicatorPtr comm);
            
            hmemo::HArray<IndexType> localIndices; //!< Global indices of the local receivers
            
            IndexType numReceiversGlobal; //!< Number of receivers global
            IndexType numReceiversLocal; //!< Number of receivers local
            
            dmemo::DistributionPtr dist_wavefield_receivers; //!< Calculated Distribution of the receivers based on the distribution of the wavefields
            dmemo::DistributionPtr no_dist_NT; //!< No distribution of the columns of the seismogram matrix
            
            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition; //!< Matrix that stores the receiver acquisition
            IndexType numParameter; //!< Number of receiver parameters given in acquisition matrix
            lama::DenseVector<IndexType> coordinates; //!< Coordinates of receivers global (1-D coordinates)
            lama::DenseVector<IndexType> receiver_type; //!< Type of Receivers: 1==Pressure, 2==vX, 3==vY, 4==vZ
            
        };
    }
}


/*! \brief Getter method for reference to receiver type
 */
template<typename ValueType>
lama::DenseVector<IndexType> const& KITGPI::Acquisition::Receivers<ValueType>::getReceiversType() const
{
    SCAI_ASSERT_ERROR( receiver_type.size() != 0 , "No receivers type set " );
    return(receiver_type);
}


/*! \brief Getter method for reference to coordinates
 */
template<typename ValueType>
lama::DenseVector<IndexType> const& KITGPI::Acquisition::Receivers<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_ERROR( coordinates.size() != 0 , "No receivers coordinates set " );
    return(coordinates);
}


/*! \brief Getter method for distribution of the receivers based on the wavefield distribution
 */
template<typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::Receivers<ValueType>::getReceiversDistribution() const
{
    SCAI_ASSERT_ERROR( dist_wavefield_receivers != nullptr , "Receivers distribution not set " );
    return(dist_wavefield_receivers);
}


/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
KITGPI::Acquisition::Receivers<ValueType>::Receivers(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
:numReceiversGlobal(0),numReceiversLocal(0),numParameter(0)
{
    init(config,ctx,dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield){
    readReceiverAcquisition(config.getReceiverFilename(),config.getNX(), config.getNY(), config.getNZ(),dist_wavefield,ctx);
}

/*! \brief Get number of global receivers
 *
 \return Number of global receivers
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversGlobal() const
{
    return(numReceiversGlobal);
}


/*! \brief Get number of local receivers
 *
 \return Number of local receivers on this process
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversLocal() const
{
    return(numReceiversLocal);
}


/*! \brief Write receivers acquisition to file
 *
 \param filename Filename to write receivers acquisition
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::writeReceiverAcquisition(std::string const& filename) const
{
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
void KITGPI::Acquisition::Receivers<ValueType>::readReceiverAcquisition(std::string const& filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx)
{

    SCAI_ASSERT_ERROR(NX>0, "NX<=0");
    SCAI_ASSERT_ERROR(NY>0, "NX<=0");
    SCAI_ASSERT_ERROR(NZ>0, "NX<=0");
    
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
        utilskernel::LArray<IndexType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<IndexType> write_coordinates_LA(*coordinates_LA);
        
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
    acquisition.getRow(receiver_type,3);
    
    /* Redistribute receiver parameter vectors to corresponding processes */
    coordinates.redistribute(dist_wavefield_receivers);
    receiver_type.redistribute(dist_wavefield_receivers);
    
    coordinates.setContextPtr(ctx);
    receiver_type.setContextPtr(ctx);
}


/*! \brief Calculation of the receiver distribution
 *
 * Calculation of the receiver distribution based on the global number of receivers numReceiversGlobal and the local number of receivers numReceiversLocal on each node.
 * Generates a GeneralDistribution.
 *
 \param comm Communicator for the generated distribution
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getReceiverDistribution(dmemo::CommunicatorPtr comm)
{
    SCAI_ASSERT(numReceiversGlobal>0," There is no global receiver (numReceiversGlobal==0)! ");
    
    dmemo::DistributionPtr dist_temp( new dmemo::GeneralDistribution(numReceiversGlobal,localIndices,comm));

    dist_wavefield_receivers=dist_temp;
}


/*! \brief Determine the number of local receivers from the global receivers coordinates
 *
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getLocalReceivers(dmemo::DistributionPtr dist_wavefield)
{
    SCAI_ASSERT(coordinates.size()>0," The vector coordinates does not contain any elements ! ");
    
    this->Global2Local(coordinates,localIndices,dist_wavefield);
    
    numReceiversLocal=localIndices.size();
    numReceiversGlobal=coordinates.size();
    
}


