#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal/all.hpp"
#include "Coordinates.hpp"
#include "Seismogram.hpp"

#include "SeismogramHandler.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        //! Handling of sources
        /*!
         * This class accounts for the handling of seismic sources.
         * It provides the reading from the source acquisition from file, the distribution of the sources and the generation of synthetic signals.
         */
        template <typename ValueType>
        class Sources
        {
            
        public:
            
            Sources():numSourcesGlobal(0),numSourcesLocal(0),numParameter(0){};
            explicit Sources(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx,dmemo::DistributionPtr dist_wavefield);
            ~Sources(){};
            
            void init(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);
            void readSourceAcquisition(std::string const& filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield);
            
            void writeSourceAcquisition(std::string const& filename) const;
            void writeSignalsToFileRaw(std::string const& filename) const;
            
            void generateSignals(IndexType NT, ValueType DT, hmemo::ContextPtr ctx);
            
            /* Getter member functions */
            lama::DenseVector<IndexType> const& getCoordinates() const;
            lama::DenseVector<IndexType> const& getSourceType() const;
            SeismogramHandler<ValueType> const& getSeismogramHandler() const;
            Seismogram<ValueType> const& getSignals() const;
            IndexType getNumSourcesGlobal() const;
            IndexType getNumSourcesLocal() const;
            
        private:
            
            dmemo::DistributionPtr getSourceDistribution(lama::DenseVector<IndexType>const& coordinates, dmemo::DistributionPtr const dist_wavefield) const;
            
            void allocateSeismogram(IndexType NT,hmemo::ContextPtr ctx);
            void generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT);
            void initSeismogramHandler(IndexType const NT,hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield);
            
            IndexType numSourcesGlobal; //!< Number of sources global
            IndexType numSourcesLocal; //!< Number of sources local
            
            dmemo::DistributionPtr dist_wavefield_sources; //!< Calculated Distribution of the sources based on the distribution of the wavefields
            
            //! Source signals
            Seismogram<ValueType> signals;
            SeismogramHandler<ValueType> sources;
            
            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition; //!< Matrix that stores the source acquisition
            IndexType numParameter; //!< Number of source parameters given in acquisition matrix
            lama::DenseVector<IndexType> coordinates; //!< Coordinates of sources global (1-D coordinates)
            lama::DenseVector<IndexType> source_type; //!< Type of source: 1==P, 2==vX, 3==vY, 4==vZ
            lama::DenseVector<IndexType> wavelet_type; //!< Type of wavelet: 1==Synthetic
            
            /* Optional acquisition Settings */
            lama::DenseVector<IndexType> wavelet_shape; //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            lama::DenseVector<ValueType> wavelet_fc; //!< Center frequency of synthetic wavelet
            lama::DenseVector<ValueType> wavelet_amp; //!< Amplitude of synthetic wavelet
            lama::DenseVector<ValueType> wavelet_tshift; //!< Time shift of synthetic wavelet
            
            
        };
    }
}

/*! \brief Get reference to SeismogramHandler
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType> const& KITGPI::Acquisition::Sources<ValueType>::getSeismogramHandler() const
{
    return(sources);
}


template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::initSeismogramHandler(IndexType const NT,hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield)
{
    IndexType const NUM_ELEMENTS_ENUM=4;
    
    IndexType count[NUM_ELEMENTS_ENUM]={0,0,0,0};
    lama::DenseVector<IndexType> coord[NUM_ELEMENTS_ENUM];
    dmemo::DistributionPtr dist[NUM_ELEMENTS_ENUM];
    
    IndexType numSourcesGlobal=source_type.size();
    
    /* Count elements for each source type */
    lama::Scalar tempScalar;
    IndexType tempIndexType;
    for(IndexType i=0; i<numSourcesGlobal;++i){
        tempScalar=source_type.getValue(i);
        tempIndexType=tempScalar.getValue<IndexType>()-1;
        
        SCAI_ASSERT_DEBUG(tempIndexType >=0 && tempIndexType <=3, "Unkown Source Type");
        ++count[tempIndexType];
    }
    
    /* Allocate lama vectors */
    for(IndexType i=0; i<NUM_ELEMENTS_ENUM; ++i){
        coord[i].allocate(count[i]);
        count[i]=0;
    }
    
    /* Sort coordinates */
    for(IndexType i=0; i<numSourcesGlobal;++i){
        
        tempScalar=source_type.getValue(i);
        tempIndexType=tempScalar.getValue<IndexType>()-1;
        
        coord[tempIndexType].setValue(count[tempIndexType], coordinates.getValue(i));
        ++count[tempIndexType];
        
    }
    
    /* Calculate distribution, redistribute coordinates and set coordinates to seismogramHandler */
    for(IndexType i=0; i<NUM_ELEMENTS_ENUM; ++i){
        if(coord[i].size()>0){
            dist[i]=getSourceDistribution(coord[i],dist_wavefield);
            sources.getSeismogram(i).allocate(ctx,dist[i],NT);
            coord[i].redistribute(dist[i]);
            sources.getSeismogram(i).setCoordinates(coord[i]);
        }
        count[i]=0;
    }
    
    lama::DenseVector<ValueType> temp;
    
    /* Copy data to the seismogram handler */
    for(IndexType i=0; i<numSourcesGlobal;++i){
        tempScalar=source_type.getValue(i);
        tempIndexType=tempScalar.getValue<IndexType>()-1;
        
        signals.getData().getRow(temp, i);
        SCAI_ASSERT_DEBUG(temp.size()==NT, "Size mismatch");
        
        sources.getSeismogram(tempIndexType).getData().setRow(temp,count[tempIndexType],utilskernel::binary::BinaryOp::COPY);
        
        ++count[tempIndexType];
    }
}


/*! \brief Get reference to source type
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
lama::DenseVector<IndexType> const& KITGPI::Acquisition::Sources<ValueType>::getSourceType() const
{
    SCAI_ASSERT_DEBUG(numSourcesGlobal == source_type.size(), "Size mismatch");
    return(source_type);
}


/*! \brief Get reference to source coordinates
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
lama::DenseVector<IndexType> const& KITGPI::Acquisition::Sources<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_DEBUG(numSourcesGlobal == coordinates.size(), "Size mismatch");
    return(coordinates);
}


/*! \brief Get reference to signals matrix
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template<typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>const& KITGPI::Acquisition::Sources<ValueType>::getSignals() const
{
    return(signals);
}


/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
KITGPI::Acquisition::Sources<ValueType>::Sources(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
:numSourcesGlobal(0),numSourcesLocal(0),numParameter(0)
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
void KITGPI::Acquisition::Sources<ValueType>::init(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
{
    readSourceAcquisition(config.getSourceFilename(),config.getNX(), config.getNY(), config.getNZ(),dist_wavefield);
    generateSignals(config.getNT(),config.getDT(),ctx);
    signals.redistribute(dist_wavefield_sources);
    initSeismogramHandler(config.getNT(),ctx,dist_wavefield);
}


/*! \brief Get number of global sources
 *
 \return Number of global sources
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Sources<ValueType>::getNumSourcesGlobal() const
{
    return(numSourcesGlobal);
}


/*! \brief Get number of local sources
 *
 \return Number of local sources on this process
 */
template<typename ValueType>
IndexType KITGPI::Acquisition::Sources<ValueType>::getNumSourcesLocal() const
{
    return(numSourcesLocal);
}


/*! \brief Write source acquisition to file
 *
 \param filename Filename to write source acquisition
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::writeSourceAcquisition(std::string const& filename) const
{
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
 \param ctx Context
 */
template<typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::readSourceAcquisition(std::string const& filename,IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield)
{
    
    SCAI_ASSERT_ERROR(NX>0, "NX<=0");
    SCAI_ASSERT_ERROR(NY>0, "NX<=0");
    SCAI_ASSERT_ERROR(NZ>0, "NX<=0");
    
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
        utilskernel::LArray<IndexType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<IndexType> write_coordinates_LA(*coordinates_LA);
        
        Coordinates<ValueType> coord;
        
        /* 2. Calculate 1-D coordinates form 3-D coordinates */
        IndexType X,Y,Z;
        for(IndexType i=0; i<numSourcesGlobal; i++){
            
            X=read_acquisition_HA[ i + numSourcesGlobal*0 ];
            Y=read_acquisition_HA[ i + numSourcesGlobal*1 ];
            Z=read_acquisition_HA[ i + numSourcesGlobal*2 ];
            
            write_coordinates_LA[i]=coord.coordinate2index(X,Y,Z,NX,NY,NZ);
        }
        
        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
        
    }
    
    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numSourcesGlobal);
    
    /* Get local sources from global sources */
    dist_wavefield_sources=getSourceDistribution(coordinates,dist_wavefield);
    
    numSourcesLocal=dist_wavefield_sources->getLocalSize();
    numSourcesGlobal=dist_wavefield_sources->getGlobalSize();
    
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
void KITGPI::Acquisition::Sources<ValueType>::writeSignalsToFileRaw(std::string const& filename) const
{
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
void KITGPI::Acquisition::Sources<ValueType>::allocateSeismogram(IndexType NT, hmemo::ContextPtr ctx)
{
    SCAI_ASSERT_DEBUG(NT>0, "NT<=0");
    if(dist_wavefield_sources==NULL) {
        COMMON_THROWEXCEPTION ( "Row distribution of sources (dist_wavefield_sources) is not set!" )
    }
    
    /* Signals matix is row distributed according to dist_wavefield_sources, No column distribution */
    signals.allocate(ctx,dist_wavefield_sources,NT);
    signals.setCoordinates(coordinates);
    signals.setTraceType(source_type);
    signals.setContextPtr(ctx);
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
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(IndexType NT, ValueType DT, hmemo::ContextPtr ctx){
    
    SCAI_ASSERT(numParameter>=5,"Number of source parameters < 5. Cannot generate signals. ");
    SCAI_ASSERT_DEBUG(NT>0, "NT<=0");
    SCAI_ASSERT_DEBUG(DT>0, "DT<=0");
    
    allocateSeismogram(NT,ctx);
    
    signals.setDT(DT);
    
    utilskernel::LArray<IndexType>* wavelet_type_LA=&wavelet_type.getLocalValues();
    hmemo::ReadAccess<IndexType> read_wavelet_type_LA(*wavelet_type_LA);
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
void KITGPI::Acquisition::Sources<ValueType>::generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT)
{
    
    SCAI_ASSERT(numParameter>=9, "Number of source parameters <= 9. Cannot generate synthetic signals. ");
    
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
    
    lama::DenseMatrix<ValueType>& signalsMatrix=signals.getData();
    signalsMatrix.setRow(signalVector,SourceLocal,utilskernel::binary::BinaryOp::COPY);
    
}


template<typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::Sources<ValueType>::getSourceDistribution(lama::DenseVector<IndexType>const& coordinates, dmemo::DistributionPtr const dist_wavefield) const
{
    SCAI_ASSERT_DEBUG(coordinates.size()>0, " The vector coordinates does not contain any elements ! ");
    
    hmemo::HArray<IndexType> localIndices;
    
    Coordinates<ValueType> coord;
    coord.Global2Local(coordinates,localIndices,dist_wavefield);
    
    dmemo::DistributionPtr dist_temp( new dmemo::GeneralDistribution(numSourcesGlobal,localIndices,dist_wavefield->getCommunicatorPtr()));
    
    return(dist_temp);
}




