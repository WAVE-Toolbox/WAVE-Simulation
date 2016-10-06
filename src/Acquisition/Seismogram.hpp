
#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Sources.hpp"
#include "Receivers.hpp"
#include "Coordinates.hpp"
#include "Configuration.hpp"

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
    void redistribute(dmemo::DistributionPtr distRow,dmemo::DistributionPtr distColumn=NULL);
    void replicate();

    void reset();
    
    //! Not yet implemented
    void normalize();
    
    IndexType getNumTraces();
    IndexType getNumTracesLocal();
    IndexType getNumSamples();
    ValueType getDT();
    
private:
    
    IndexType numSamples=0; //!< Number of samples of one trace
    IndexType numTracesGlobal=0; //!< Number of global traces
    IndexType numTracesLocal=0; //!< Number of local traces
    
    /* header information */
    ValueType DT=0; //!< Temporal sampling in seconds
    ValueType receiver_type; //!< Type of receiver
    
    /* raw data */
    lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
    
};

//! \brief Replicate seismogram on all processes
/*!
 * Creates a copy of the seismogram on all processe
 */
template <typename ValueType>
void Seismogram<ValueType>::replicate()
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
 \param distSeismogram Distribution for traces
 \param NT Total number of samples per trace
 */
template <typename ValueType>
void Seismogram<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT)
{
    dmemo::DistributionPtr no_dist_NT( new scai::dmemo::NoDistribution ( NT ) );
    
    data.setContext(ctx);
    
    data.allocate(distSeismogram,no_dist_NT);
}


//! \brief reset seismogram set the seismogram data to zero
/*!
 */
template <typename ValueType>
void Seismogram<ValueType>::reset()
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
void Seismogram<ValueType>::redistribute(dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
{
    if(distColumn==NULL){
        dmemo::DistributionPtr distSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    }
    
    data.redistribute(distTraces,distSamples);
}


//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void Seismogram<ValueType>::ReadFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
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
void Seismogram<ValueType>::writeToFileRaw(std::string filename)
{
    data.writeToFile(filename);
}


//! \brief Get temporal sampling
template <typename ValueType>
ValueType Seismogram<ValueType>::getDT()
{
    return(DT);
}


//! \brief Get number of samples per trace
template <typename ValueType>
IndexType Seismogram<ValueType>::getNumSamples()
{
    return(numSamples);
}


//! \brief Get number of local traces
template <typename ValueType>
IndexType Seismogram<ValueType>::getNumTracesLocal()
{
    return(numTracesLocal);
}

//! \brief Get number of global traces
template <typename ValueType>
IndexType Seismogram<ValueType>::getNumTraces()
{
    return(numTracesGlobal);
}



