
#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>
#include "Sourcesignal.hpp"
#include "Coordinates.hpp"
#include "Configuration.hpp"

#pragma once

template <typename ValueType>
class Sources : private Sourcesignal<ValueType>, private Coordinates<ValueType>
{
    
public:
    
    Sources(){};
    Sources(Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield);
    ~Sources(){};
    
    void readSourceAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield);
    void generateSignals(IndexType NT, ValueType DT);
    void applySource(lama::DenseVector<ValueType>& wavefield, IndexType nt);

    void writeSignals(std::string filename);
    void writeSourceAcquisition(std::string filename);
    
    IndexType getNumSourcesGlobal();
    IndexType getNumSourcesLocal();
    
private:

    void getLocalSources(dmemo::DistributionPtr dist_wavefield);
    void getSourceDistribution(dmemo::CommunicatorPtr comm);
    void allocateSignals(IndexType NT);
    void generateSyntheticSignal(IndexType nsourcelocal_this, IndexType NT, ValueType DT);
    
    IndexType nsourcesglobal=0; //!< Number of sources global
    IndexType nsourceslocal=0; //!< Number of sources local
    
    lama::DenseVector<ValueType> coordinates; //!< Coordinates of sources global (1-D coordinates)

    dmemo::DistributionPtr dist_wavefield_sources=NULL;
    dmemo::DistributionPtr no_dist_NT=NULL;
    
    //! Source signals
    lama::DenseMatrix<ValueType> signals;
    
    /* Acquisition Settings */
    lama::DenseMatrix<ValueType> acquisition; //!< Matrix that stores the source acquisition
    IndexType nparameter=0; //!< Number of source parameters given in acquisition matrix
    lama::DenseVector<ValueType> source_type; //!< Type of source: 1==Pressure
    lama::DenseVector<ValueType> wavelet_type; //!< Type of wavelet: 1==Synthetic
    
    /* Optional acquisition Settings */
    lama::DenseVector<ValueType> wavelet_shape; //!< Shape of wavelet: 1==Ricker
    lama::DenseVector<ValueType> wavelet_fc; //!< Center frequency of synthetic wavelet
    lama::DenseVector<ValueType> wavelet_amp; //!< Amplitude of synthetic wavelet
    lama::DenseVector<ValueType> wavelet_tshift; //!< Time shift of synthetic wavelet
    
    
};


template<typename ValueType>
Sources<ValueType>::Sources(Configuration<ValueType> config, dmemo::DistributionPtr dist_wavefield){
    readSourceAcquisition(config.getSourceFilename(),config.getNX(), config.getNY(), config.getNZ(),dist_wavefield);
    generateSignals(config.getNT(),config.getDT());
}

template<typename ValueType>
IndexType Sources<ValueType>::getNumSourcesGlobal(){
    return(nsourcesglobal);
}

template<typename ValueType>
IndexType Sources<ValueType>::getNumSourcesLocal(){
    return(nsourceslocal);
}

template<typename ValueType>
void Sources<ValueType>::writeSourceAcquisition(std::string filename){
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.assignTranspose(acquisition);
    acquisition_temp.writeToFile(filename);
}

template<typename ValueType>
void Sources<ValueType>::readSourceAcquisition(std::string filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield)
{
    
    /* Read acquisition matrix */
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(filename);
    
    IndexType nrow_temp=acquisition_temp.getNumRows();
    IndexType ncolumn_temp=acquisition_temp.getNumColumns();
    
    /* Derive number of sources and number of read-in parameters */
    nsourcesglobal=nrow_temp;
    nparameter=ncolumn_temp;
    
    /* Check if number of parameters is supported */
    if( nparameter<5 || nparameter>9 ) {
        COMMON_THROWEXCEPTION ( "Source acquisition file has an unkown format " )
    }
    
    /* Distribution: Master process only (comm->myRank()==0) */
    dmemo::DistributionPtr dist_master_nparameter( new dmemo::CyclicDistribution( nparameter, nparameter, dist_wavefield->getCommunicatorPtr() ) );
    dmemo::DistributionPtr dist_master_nsourcesglobal( new dmemo::CyclicDistribution( nsourcesglobal, nsourcesglobal, dist_wavefield->getCommunicatorPtr() )  );

    /* Distribution: Replicated on all processes */
    dmemo::DistributionPtr no_dist_nsourcesglobal( new scai::dmemo::NoDistribution ( nsourcesglobal ) );
    dmemo::DistributionPtr no_dist_nparameter( new scai::dmemo::NoDistribution ( nparameter ) );

    /* Allocate acquisition matrix on master */
    acquisition.allocate(dist_master_nparameter,no_dist_nsourcesglobal);
    
    /* Allocate coordinates on master */
    coordinates.allocate(dist_master_nsourcesglobal);
    
    /* Allocate source parameter vectors on master */
    source_type.allocate(dist_master_nsourcesglobal);
    wavelet_type.allocate(dist_master_nsourcesglobal);
    if(nparameter>5){
        wavelet_shape.allocate(dist_master_nsourcesglobal);
        wavelet_fc.allocate(dist_master_nsourcesglobal);
        wavelet_amp.allocate(dist_master_nsourcesglobal);
        wavelet_tshift.allocate(dist_master_nsourcesglobal);
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
        for(IndexType i=0; i<nsourcesglobal; i++){
            
            X=read_acquisition_HA[ i + nsourcesglobal*0 ];
            Y=read_acquisition_HA[ i + nsourcesglobal*1 ];
            Z=read_acquisition_HA[ i + nsourcesglobal*2 ];
    
            write_coordinates_LA[i]=this->coordinate2index(X,Y,Z,NX,NY,NZ);
        }
        
        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
        
    }
    
    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_nsourcesglobal);
    
    /* Get local sources from global sources */
    getLocalSources(dist_wavefield);
    
    /* Get source distribution */
    getSourceDistribution(dist_wavefield->getCommunicatorPtr());

    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_nparameter,no_dist_nsourcesglobal);
    
    /* Allocate source parameter vectors on all processes */
    source_type.allocate(nsourcesglobal);
    wavelet_type.allocate(nsourcesglobal);
    if(nparameter>5){
        wavelet_shape.allocate(nsourcesglobal);
        wavelet_fc.allocate(nsourcesglobal);
        wavelet_amp.allocate(nsourcesglobal);
        wavelet_tshift.allocate(nsourcesglobal);
    }
    
    /* Save source configurations from acquisition matrix in vectors */
    acquisition.getLocalRow(source_type,3);
    acquisition.getLocalRow(wavelet_type,4);
    if(nparameter>5){
        acquisition.getLocalRow(wavelet_shape,5);
        acquisition.getLocalRow(wavelet_fc,6);
        acquisition.getLocalRow(wavelet_amp,7);
        acquisition.getLocalRow(wavelet_tshift,8);
    }

    /* Redistribute source parameter vectors to corresponding processes */
    coordinates.redistribute(dist_wavefield_sources);
    source_type.redistribute(dist_wavefield_sources);
    wavelet_type.redistribute(dist_wavefield_sources);
    if(nparameter>5){
        wavelet_shape.redistribute(dist_wavefield_sources);
        wavelet_fc.redistribute(dist_wavefield_sources);
        wavelet_amp.redistribute(dist_wavefield_sources);
        wavelet_tshift.redistribute(dist_wavefield_sources);
    }
    
}


template<typename ValueType>
void Sources<ValueType>::writeSignals(std::string filename){
    signals.writeToFile(filename);
}

template<typename ValueType>
void Sources<ValueType>::allocateSignals(IndexType NT)
{
    dmemo::DistributionPtr no_dist_NT_temp( new scai::dmemo::NoDistribution ( NT ) );
    no_dist_NT=no_dist_NT_temp;
    
    if(dist_wavefield_sources==NULL) {
        COMMON_THROWEXCEPTION ( "Row distribution of sources is not set!" )
    }
    
    /* Signals matix is row distributed according to dist_wavefield_sources, No column distribution */
    signals.allocate(dist_wavefield_sources,no_dist_NT);
}


template <typename ValueType>
void Sources<ValueType>::generateSignals(IndexType NT, ValueType DT){
    
    allocateSignals(NT);
    
    utilskernel::LArray<ValueType>* wavelet_type_LA=&wavelet_type.getLocalValues();
    hmemo::ReadAccess<ValueType> read_wavelet_type_LA(*wavelet_type_LA);
    IndexType wavelet_type_i;
    
    for(IndexType i=0; i<nsourceslocal; i++){
        
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


template <typename ValueType>
void Sources<ValueType>::generateSyntheticSignal(IndexType nsourcelocal_this, IndexType NT, ValueType DT){
    
    lama::DenseVector<ValueType> signal;
    signal.allocate(NT);
    
    /* Cast to IndexType */
    IndexType wavelet_shape_i=wavelet_shape.getLocalValues()[nsourcelocal_this];
    
    switch (wavelet_shape_i) {
        case 1:
            /* Ricker */
            this->Ricker(signal,  NT,  DT,  wavelet_fc.getLocalValues()[nsourcelocal_this],  wavelet_amp.getLocalValues()[nsourcelocal_this],  wavelet_tshift.getLocalValues()[nsourcelocal_this]);
            break;
            
        default:
             COMMON_THROWEXCEPTION ( "Unkown wavelet shape ")
            break;
    }
    
    utilskernel::LArray<ValueType>* signal_LA=&signal.getLocalValues();
    hmemo::ReadAccess<ValueType> read_signal(*signal_LA);
    
    lama::DenseStorage<ValueType>* signals_DS=&signals.getLocalStorage();
    hmemo::HArray<ValueType>* signals_HA=&signals_DS->getData();
    hmemo::WriteAccess<ValueType> write_signals_HA(*signals_HA);
    
    for(IndexType i=0; i<NT; i++){
        write_signals_HA[i+NT*nsourcelocal_this]=read_signal[i];
    }
    read_signal.release();
    write_signals_HA.release();
    
}


template<typename ValueType>
void Sources<ValueType>::getSourceDistribution(dmemo::CommunicatorPtr comm)
{
    dmemo::DistributionPtr dist_temp( new dmemo::GenBlockDistribution(nsourcesglobal,nsourceslocal,comm));
    dist_wavefield_sources=dist_temp;
}



/*! \brief Determine local source coordinates from global distribution
 * Determines the local source coordinates from the global distribution
 *
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void Sources<ValueType>::getLocalSources(dmemo::DistributionPtr dist_wavefield)
{
    lama::DenseVector<ValueType> coordinateslocal; //!< Coordinates of sources local (1-D coordinates)
    
    this->Global2Local(coordinates,coordinateslocal,dist_wavefield);
    
    nsourceslocal=coordinateslocal.size();
    nsourcesglobal=coordinates.size();
    
}


template<typename ValueType>
void Sources<ValueType>::applySource(lama::DenseVector<ValueType>& wavefield, IndexType nt)
{
    scai::lama::Scalar source_index_temp;
    IndexType source_index;
    
    for(IndexType i=0; i<nsourcesglobal; i++){
        
        source_index_temp=coordinates.getValue(i);
        source_index=source_index_temp.getValue<IndexType>();
        
        wavefield.setValue(source_index, wavefield.getValue(source_index) + signals.getValue(i,nt) );
        
    }
    
}


