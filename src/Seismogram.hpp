
#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Sources.hpp"
#include "Receivers.hpp"
#include "Coordinates.hpp"
#include "Configuration.hpp"


template <typename ValueType>
class Seismogram
{
    
public:
    
    Seismogram(){};
    
    ~Seismogram(){};
    
    void writeToFileRaw(std::string filename);
    void ReadFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples);
    
    void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT);
    void redistribute(dmemo::DistributionPtr distRow,dmemo::DistributionPtr distColumn=NULL);
    void replicate();

    void reset();
    
    void normalize();
    
    IndexType getNumTraces();
    IndexType getNumTracesLocal();
    IndexType getNumSamples();
    ValueType getDT();
    
private:
    
    IndexType numSamples=0;
    IndexType numTracesGlobal=0;
    IndexType numTracesLocal=0;
    
    /* header information */
    ValueType DT=0;
    ValueType receiver_type;
    lama::DenseVector<ValueType> offset;
    
    /* raw data */
    lama::DenseMatrix<ValueType> data;
    
};


template <typename ValueType>
void Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    
    dmemo::DistributionPtr no_dist_Traces( new scai::dmemo::NoDistribution ( numTracesGlobal ) );
    
    redistribute(no_dist_Traces,no_dist_numSamples);
}

template <typename ValueType>
void Seismogram<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT)
{
    dmemo::DistributionPtr no_dist_NT( new scai::dmemo::NoDistribution ( NT ) );
    
    data.setContext(ctx);
    
    data.allocate(distSeismogram,no_dist_NT);
}


template <typename ValueType>
void Seismogram<ValueType>::reset()
{
    data.scale(0.0);
}

template <typename ValueType>
void Seismogram<ValueType>::redistribute(dmemo::DistributionPtr distRow,dmemo::DistributionPtr distColumn)
{
    if(distColumn==NULL){
        dmemo::DistributionPtr distColumn( new scai::dmemo::NoDistribution ( numSamples ) );
    }
    
    data.redistribute(distRow,distColumn);
}


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

template <typename ValueType>
void Seismogram<ValueType>::writeToFileRaw(std::string filename)
{
    data.writeToFile(filename);
}

template <typename ValueType>
ValueType Seismogram<ValueType>::getDT()
{
    return(DT);
}

template <typename ValueType>
IndexType Seismogram<ValueType>::getNumSamples()
{
    return(numSamples);
}

template <typename ValueType>
IndexType Seismogram<ValueType>::getNumTracesLocal()
{
    return(numTracesLocal);
}


template <typename ValueType>
IndexType Seismogram<ValueType>::getNumTraces()
{
    return(numTracesGlobal);
}



