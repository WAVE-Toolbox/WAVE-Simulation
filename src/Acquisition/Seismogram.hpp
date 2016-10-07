
#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Sources.hpp"
#include "Receivers.hpp"
#include "Coordinates.hpp"
#include "../Configuration/Configuration.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        //! Seismogram class
        /*!
         * This class handels a seismogram which consists of several traces.
         */
        template <typename ValueType>
        class Seismogram
        {
            
        public:
            
            //! Default constructor
            Seismogram(){};
            
            //! Default destructor
            ~Seismogram(){};
            
            void writeToFileRaw(std::string filename);
            void ReadFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples);
            
            void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT);
            
            void init(Receivers<ValueType>& receiver, IndexType NT, hmemo::ContextPtr ctx);
            
            void redistribute(dmemo::DistributionPtr distRow,dmemo::DistributionPtr distColumn=NULL);
            void replicate();
            
            void resetData();
            
            //! Not yet implemented
            void normalize();
            
            IndexType getNumTraces();
            IndexType getNumTracesLocal();
            IndexType getNumSamples();
            ValueType getDT();
            
            lama::DenseMatrix<ValueType>* getData();
            lama::DenseVector<ValueType>* getReceiverType();
            lama::DenseVector<ValueType>* getCoordinates();
            
        private:
            
            IndexType numSamples=0; //!< Number of samples of one trace
            IndexType numTracesGlobal=0; //!< Number of global traces
            IndexType numTracesLocal=0; //!< Number of local traces
            
            /* header information */
            ValueType DT=0.0; //!< Temporal sampling in seconds
            lama::DenseVector<ValueType> receiver_type; //!< Type of receiver
            lama::DenseVector<ValueType> coordinates; //!< Coordinates of the traces
            
            /* raw data */
            lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
            
        };
    }
}


//! \brief Initiate the seismogram by receivers
/*!
 *
 \param receiver Receivers which will record into this seismogram
 \param NT Total number of time steps
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::init(Receivers<ValueType>& receiver, IndexType NT, hmemo::ContextPtr ctx){
    
    /* Allocation */
    dmemo::DistributionPtr dist_traces=receiver.getReceiversDistribution();
    allocate(ctx,dist_traces,NT);
    
    /* set header information */
    lama::DenseVector<ValueType>& coordinates_temp=*receiver.getCoordinates();
    lama::DenseVector<ValueType>& receiver_type_temp=*receiver.getReceiversType();
    coordinates=coordinates_temp;
    receiver_type=receiver_type_temp;
    
    numTracesLocal=receiver.getNumReceiversLocal();
    numTracesGlobal=receiver.getNumReceiversGlobal();
    numSamples=NT;
}


//! \brief Get reference to Receiver Type
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<ValueType>* KITGPI::Acquisition::Seismogram<ValueType>::getReceiverType(){
    return(&receiver_type);
}


//! \brief Get reference to coordinates
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<ValueType>* KITGPI::Acquisition::Seismogram<ValueType>::getCoordinates(){
    return(&coordinates);
}


//! \brief Get reference to seismogram data
/*!
 * 
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType>* KITGPI::Acquisition::Seismogram<ValueType>::getData(){
    return(&data);
}


//! \brief Replicate seismogram on all processes
/*!
 * Creates a copy of the seismogram on all processe
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    
    dmemo::DistributionPtr no_dist_Traces( new scai::dmemo::NoDistribution ( numTracesGlobal ) );
    
    redistribute(no_dist_Traces,no_dist_numSamples);
}


//! \brief Allocate seismogram
/*!
 * Allocates seismogram based on a given distribution of the traces and the number of samples per trace
 *
 \param ctx Context
 \param distTraces Distribution for traces
 \param NT Total number of samples per trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distTraces, IndexType NT)
{
    
    numSamples=NT;
    numTracesGlobal=distTraces->getGlobalSize();
    numTracesLocal=distTraces->getLocalSize();
    
    dmemo::DistributionPtr no_dist_NT( new scai::dmemo::NoDistribution ( NT ) );
    
    data.setContextPtr(ctx);
    receiver_type.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);
    
    data.allocate(distTraces,no_dist_NT);
    receiver_type.allocate(distTraces);
    coordinates.allocate(distTraces);
}


//! \brief reset seismogram set the seismogram data to zero
/*!
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetData()
{
    data.scale(0.0);
}


//! \brief Redistribute seismogram data
/*!
 *
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::redistribute(dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
{
    if(distSamples==NULL){
        dmemo::DistributionPtr distSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    }
    
    data.redistribute(distTraces,distSamples);
    receiver_type.redistribute(distTraces);
    coordinates.redistribute(distTraces);
}


//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::ReadFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
{
    data.ReadFromFile(filename);
    
    IndexType nrow_temp=data.getNumRows();
    IndexType ncolumn_temp=data.getNumColumns();
    
    numSamples=ncolumn_temp;
    numTracesGlobal=nrow_temp;
    
    if(distTraces==NULL && distSamples==NULL){
        replicate();
    } else {
        redistribute(distTraces,distSamples);
    }
    
}


//! \brief Write a seismogram to disk without header
/*!
 *
 \param filename Filename to write seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileRaw(std::string filename)
{
    if(data.getNumValues()==0) {
        COMMON_THROWEXCEPTION("Seismogramm data is not allocated")
    }
    data.writeToFile(filename);
}


//! \brief Get temporal sampling
template <typename ValueType>
ValueType KITGPI::Acquisition::Seismogram<ValueType>::getDT()
{
    if(DT==0.0){
        COMMON_THROWEXCEPTION("DT is not set for this seismogram")
    }
    return(DT);
}


//! \brief Get number of samples per trace
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumSamples()
{
    if(numSamples==0){
        numSamples=data.getNumColumns();
        if(numSamples==0){
            COMMON_THROWEXCEPTION("Seismogram is not allocated")
        }
    }
    return(numSamples);
}


//! \brief Get number of local traces
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesLocal()
{
    return(numTracesLocal);
}


//! \brief Get number of global traces
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTraces()
{
    if(numTracesGlobal==0){
        numTracesGlobal=data.getNumRows();
        if(numTracesGlobal==0){
            COMMON_THROWEXCEPTION("Seismogram is not allocated")
        }
    }
    return(numTracesGlobal);
}



