#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal/all.hpp"
#include "Coordinates.hpp"
#include "Seismogram.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        //! Handling of sources
        /*!
         * This class accounts for the handling of seismic sources.
         * It provides the reading from the source acquisition from file, the distribution of the sources and the generation of synthetic signals.
         */
        template <typename ValueType>
        class Sources : private Coordinates<ValueType>
        {
            
        public:
            
            Sources():numSourcesGlobal(0),numSourcesLocal(0),numParameter(0){};
            Sources(Configuration::Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield);
            ~Sources(){};
            
            void readSourceAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield);
            void writeSourceAcquisition(std::string filename);
            
            void generateSignals(IndexType NT, ValueType DT);
            void writeSignalsToFileRaw(std::string filename);
            
            lama::DenseVector<ValueType>& getCoordinates();
            lama::DenseVector<ValueType>& getSourceType();
            Seismogram<ValueType>& getSignals();
            
            IndexType getNumSourcesGlobal();
            IndexType getNumSourcesLocal();
            
        private:
            
            void getLocalSources(dmemo::DistributionPtr dist_wavefield);
            void getSourceDistribution(dmemo::CommunicatorPtr comm);
            void allocateSignals(IndexType NT);
            void generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT);
            
            IndexType numSourcesGlobal; //!< Number of sources global
            IndexType numSourcesLocal; //!< Number of sources local
            
            dmemo::DistributionPtr dist_wavefield_sources; //!< Calculated Distribution of the sources based on the distribution of the wavefields
            dmemo::DistributionPtr no_dist_NT; //!< No distribution of the columns of the signals matrix
            
            hmemo::HArray<IndexType> localIndices; //!< Global indices of the local sources
            
            //! Source signals
            Seismogram<ValueType> signals;
            
            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition; //!< Matrix that stores the source acquisition
            IndexType numParameter; //!< Number of source parameters given in acquisition matrix
            lama::DenseVector<ValueType> coordinates; //!< Coordinates of sources global (1-D coordinates)
            lama::DenseVector<ValueType> source_type; //!< Type of source: 1==P, 2==vX, 3==vY, 4==vZ
            lama::DenseVector<ValueType> wavelet_type; //!< Type of wavelet: 1==Synthetic
            
            /* Optional acquisition Settings */
            lama::DenseVector<ValueType> wavelet_shape; //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            lama::DenseVector<ValueType> wavelet_fc; //!< Center frequency of synthetic wavelet
            lama::DenseVector<ValueType> wavelet_amp; //!< Amplitude of synthetic wavelet
            lama::DenseVector<ValueType> wavelet_tshift; //!< Time shift of synthetic wavelet
            
        };
    }
}


/*! \brief Get reference to source type
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Acquisition::Sources<ValueType>::getSourceType(){
    return(source_type);
}


/*! \brief Get reference to source coordinates
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Acquisition::Sources<ValueType>::getCoordinates(){
    return(coordinates);
}


/*! \brief Get reference to signals matrix
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>& KITGPI::Acquisition::Sources<ValueType>::getSignals(){
    return(signals);
}


/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
KITGPI::Acquisition::Sources<ValueType>::Sources(Configuration::Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield)
:numSourcesGlobal(0),numSourcesLocal(0),numParameter(0)
{
    readSourceAcquisition(config.getSourceFilename(),config.getNX(), config.getNY(), config.getNZ(),dist_wavefield);
    generateSignals(config.getNT(),config.getDT());
}


/*! \brief Get number of global sources
 *
 \return Number of global sources
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Sources<ValueType>::getNumSourcesGlobal(){
    if(numSourcesGlobal==0){
        numSourcesGlobal=signals.getNumTracesGlobal();
        if(numSourcesGlobal==0){
            COMMON_THROWEXCEPTION("The signals are not allocated")
        }
    }
    return(numSourcesGlobal);
}


/*! \brief Get number of local sources
 *
 \return Number of local sources on this process
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Sources<ValueType>::getNumSourcesLocal(){
    return(numSourcesLocal);
}


/*! \brief Write source acquisition to file
 *
 \param filename Filename to write source acquisition
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::writeSourceAcquisition(std::string filename){
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.assignTranspose(acquisition);
    acquisition_temp.writeToFile(filename);
}


/*! \brief Read source acquisition from file
 *
 * This method reads in the source acquisition, calculates the 1-D coordinates from the 3-D coordinates,
 * splits up the source configuration into the corresponding vectors and calculates the source distribution.
 * The parameter vectors will be distributed accordingly to the source distribution.
 *
 \param filename Filename to read source acquisition
 \param NX Number of global grid points in X
 \param NY Number of global grid points in Y
 \param NZ Number of global grid points in Z
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::readSourceAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield)
{
    
    /* Read acquisition matrix */
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(filename);
    
    IndexType nrow_temp=acquisition_temp.getNumRows();
    IndexType ncolumn_temp=acquisition_temp.getNumColumns();
    
    /* Derive number of sources and number of read-in parameters */
    numSourcesGlobal=nrow_temp;
    numParameter=ncolumn_temp;
    
    /* Check if number of parameters is supported */
    if( numParameter<5 || numParameter>9 ) {
        COMMON_THROWEXCEPTION ( "Source acquisition file has an unkown format " )
    }
    
    /* Distribution: Master process only (comm->myRank()==0) */
    dmemo::DistributionPtr dist_master_numParameter( new dmemo::CyclicDistribution( numParameter, numParameter, dist_wavefield->getCommunicatorPtr() ) );
    dmemo::DistributionPtr dist_master_numSourcesGlobal( new dmemo::CyclicDistribution( numSourcesGlobal, numSourcesGlobal, dist_wavefield->getCommunicatorPtr() )  );
    
    /* Distribution: Replicated on all processes */
    dmemo::DistributionPtr no_dist_numSourcesGlobal( new scai::dmemo::NoDistribution ( numSourcesGlobal ) );
    dmemo::DistributionPtr no_dist_numParameter( new scai::dmemo::NoDistribution ( numParameter ) );
    
    /* Allocate acquisition matrix on master */
    acquisition.allocate(dist_master_numParameter,no_dist_numSourcesGlobal);
    
    /* Allocate coordinates on master */
    coordinates.allocate(dist_master_numSourcesGlobal);
    
    /* Allocate source parameter vectors on master */
    source_type.allocate(dist_master_numSourcesGlobal);
    wavelet_type.allocate(dist_master_numSourcesGlobal);
    if(numParameter>5){
        wavelet_shape.allocate(dist_master_numSourcesGlobal);
        wavelet_fc.allocate(dist_master_numSourcesGlobal);
        wavelet_amp.allocate(dist_master_numSourcesGlobal);
        wavelet_tshift.allocate(dist_master_numSourcesGlobal);
    }
    
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
        for(IndexType i=0; i<numSourcesGlobal; i++){
            
            X=read_acquisition_HA[ i + numSourcesGlobal*0 ];
            Y=read_acquisition_HA[ i + numSourcesGlobal*1 ];
            Z=read_acquisition_HA[ i + numSourcesGlobal*2 ];
            
            write_coordinates_LA[i]=this->coordinate2index(X,Y,Z,NX,NY,NZ);
        }
        
        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
        
    }
    
    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numSourcesGlobal);
    
    /* Get local sources from global sources */
    getLocalSources(dist_wavefield);
    
    /* Get source distribution */
    getSourceDistribution(dist_wavefield->getCommunicatorPtr());
    
    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_numParameter,no_dist_numSourcesGlobal);
    
    /* Allocate source parameter vectors on all processes */
    source_type.allocate(numSourcesGlobal);
    wavelet_type.allocate(numSourcesGlobal);
    if(numParameter>5){
        wavelet_shape.allocate(numSourcesGlobal);
        wavelet_fc.allocate(numSourcesGlobal);
        wavelet_amp.allocate(numSourcesGlobal);
        wavelet_tshift.allocate(numSourcesGlobal);
    }
    
    /* Save source configurations from acquisition matrix in vectors */
    acquisition.getRow(source_type,3);
    acquisition.getRow(wavelet_type,4);
    if(numParameter>5){
        acquisition.getRow(wavelet_shape,5);
        acquisition.getRow(wavelet_fc,6);
        acquisition.getRow(wavelet_amp,7);
        acquisition.getRow(wavelet_tshift,8);
    }
    
    /* Redistribute source parameter vectors to corresponding processes */
    coordinates.redistribute(dist_wavefield_sources);
    source_type.redistribute(dist_wavefield_sources);
    wavelet_type.redistribute(dist_wavefield_sources);
    if(numParameter>5){
        wavelet_shape.redistribute(dist_wavefield_sources);
        wavelet_fc.redistribute(dist_wavefield_sources);
        wavelet_amp.redistribute(dist_wavefield_sources);
        wavelet_tshift.redistribute(dist_wavefield_sources);
    }
    
}


/*! \brief Write source signals to file
 *
 \param filename Filename to write source signals
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::writeSignalsToFileRaw(std::string filename){
    signals.writeToFileRaw(filename);
}


/*! \brief Allocation of the source signals matrix
 *
 * Allocation of the source signals matrix based on an already defined source distribution and the number of time steps.
 * The source signal matrix is allocated based on the distributions.
 *
 \param NT Number of time steps
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::allocateSignals(IndexType NT)
{
    if(dist_wavefield_sources==NULL) {
        COMMON_THROWEXCEPTION ( "Row distribution of sources (dist_wavefield_sources) is not set!" )
    }
    
    /* Signals matix is row distributed according to dist_wavefield_sources, No column distribution */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();
    signals.allocate(ctx,dist_wavefield_sources,NT);
    
    signals.getCoordinates()=coordinates;
    signals.getTraceType()=source_type;
}


/*! \brief Generation of the source signals
 *
 * Allocation and calculation of the source signals accordingly to the source parameter vectors.
 * The calculation is performed locally on each node.
 *
 \param NT Number of time steps
 \param DT Time step interval
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(IndexType NT, ValueType DT){
    
    if(numParameter<5) {
        COMMON_THROWEXCEPTION ( "Number of source parameters < 5. Cannot generate signals. " )
    }
    
    allocateSignals(NT);
    
    utilskernel::LArray<ValueType>* wavelet_type_LA=&wavelet_type.getLocalValues();
    hmemo::ReadAccess<ValueType> read_wavelet_type_LA(*wavelet_type_LA);
    IndexType wavelet_type_i;
    
    for(IndexType i=0; i<numSourcesLocal; i++){
        
        /* Cast to IndexType */
        wavelet_type_i=read_wavelet_type_LA[i];
        
        switch (wavelet_type_i) {
            case 1:
                /* Synthetic wavelet */
                generateSyntheticSignal(i,NT,DT);
                break;
                
            default:
                COMMON_THROWEXCEPTION ( "Unkown wavelet type ")
                break;
        }
        
    }
    
}


/*! \brief Generation of synthetic source signals
 *
 * Calculation of a synthetic source signal accordingly to the source parameter vectors for the given local source number.
 * Uses the entries of the wavelet_shape vector to determine the shape of the wavelet.
 *
 \param SourceLocal Number of the local source
 \param NT Number of time steps
 \param DT Time step interval
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT){
    
    if(numParameter<9) {
        COMMON_THROWEXCEPTION ( "Number of source parameters <= 9. Cannot generate synthetic signals. " )
    }
    
    lama::DenseVector<ValueType> signalVector;
    signalVector.allocate(NT);
    
    /* Cast to IndexType */
    IndexType wavelet_shape_i=wavelet_shape.getLocalValues()[SourceLocal];
    
    switch (wavelet_shape_i) {
        case 1:
            /* Ricker */
            SourceSignal::Ricker<ValueType>(signalVector,  NT,  DT,  wavelet_fc.getLocalValues()[SourceLocal],  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        case 2:
            /* combination of sin signals */
            SourceSignal::SinW<ValueType>(signalVector,  NT,  DT,  wavelet_fc.getLocalValues()[SourceLocal],  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        case 3:
            /* sin3 signal */
            SourceSignal::SinThree<ValueType>(signalVector,  NT,  DT,  wavelet_fc.getLocalValues()[SourceLocal],  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        case 4:
            /* First derivative of a Gaussian (FGaussian) */
            SourceSignal::FGaussian<ValueType>(signalVector,  NT,  DT,  wavelet_fc.getLocalValues()[SourceLocal],  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        case 5:
            /* Spike signal */
            SourceSignal::Spike<ValueType>(signalVector,  NT,  DT, 0,  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        case 6:
            /* integral sin3 signal */
            SourceSignal::IntgSinThree<ValueType>(signalVector,  NT,  DT,  wavelet_fc.getLocalValues()[SourceLocal],  wavelet_amp.getLocalValues()[SourceLocal],  wavelet_tshift.getLocalValues()[SourceLocal]);
            break;
            
        default:
            COMMON_THROWEXCEPTION ( "Unkown wavelet shape ")
            break;
    }
    
    utilskernel::LArray<ValueType>* signalVector_LA=&signalVector.getLocalValues();
    hmemo::ReadAccess<ValueType> read_signalVector(*signalVector_LA);
    
    lama::DenseMatrix<ValueType>& signalsMatrix=signals.getData();
    lama::DenseStorage<ValueType>* signalsMatrix_DS=&signalsMatrix.getLocalStorage();
    hmemo::HArray<ValueType>* signalsMatrix_HA=&signalsMatrix_DS->getData();
    hmemo::WriteAccess<ValueType> write_signalsMatrix_HA(*signalsMatrix_HA);
    
    for(IndexType i=0; i<NT; i++){
        write_signalsMatrix_HA[i+NT*SourceLocal]=read_signalVector[i];
    }
    read_signalVector.release();
    write_signalsMatrix_HA.release();
    
}


/*! \brief Calculation of the source distribution
 *
 * Calculation of the source distribution based on the global number of sources numSourcesGlobal and the local number of sources numSourcesLocal on each node.
 * Generates a GeneralDistribution.
 *
 \param comm Communicator for the generated distribution
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::getSourceDistribution(dmemo::CommunicatorPtr comm)
{
    if(numSourcesGlobal==0){
        COMMON_THROWEXCEPTION ( " There is no global source (numSourcesGlobal==0)! ")
    }
    
    //dmemo::DistributionPtr dist_temp( new dmemo::GenBlockDistribution(numSourcesGlobal,numSourcesLocal,comm));
    dmemo::DistributionPtr dist_temp( new dmemo::GeneralDistribution(numSourcesGlobal,localIndices,comm));

    dist_wavefield_sources=dist_temp;
}


/*! \brief Determine the number of local source from the global source coordinates
 *
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::getLocalSources(dmemo::DistributionPtr dist_wavefield)
{
    
    if(coordinates.size()==0){
        COMMON_THROWEXCEPTION ( " The vector coordinates does not contain any elements ! ")
    }
    
    
    this->Global2Local(coordinates,localIndices,dist_wavefield);
    
    numSourcesLocal=localIndices.size();
    numSourcesGlobal=coordinates.size();
    
}



