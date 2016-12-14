#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"



namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief Derivatives namespace
        namespace SourceReceiverImpl {
            
            template<typename ValueType>
            class SourceReceiverImpl
            {
            public:
                
                SourceReceiverImpl()=delete;                
                explicit SourceReceiverImpl(Acquisition::Sources<ValueType> const& sourcesIN,Acquisition::Receivers<ValueType>& receiversIN, Wavefields::Wavefields<ValueType>& wavefieldIN);
                ~SourceReceiverImpl(){};
                
                void applySource(IndexType t);
                void gatherSeismogram(IndexType t);
                
            protected:
                
                virtual void applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) = 0;
                inline virtual void applySourceVX(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;
                inline virtual void applySourceVY(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;
                inline virtual void applySourceVZ(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;
                inline void applySourceSingle(Acquisition::Seismogram<ValueType> const& seismo, lama::DenseVector<ValueType>& wavefieldSingle, lama::DenseVector<ValueType>& temp, IndexType t);
                
                virtual void gatherSeismogramPressure(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)=0;
                inline virtual void gatherSeismogramVX(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t);
                inline virtual void gatherSeismogramVY(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t);
                inline virtual void gatherSeismogramVZ(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t);
                inline void gatherSeismogramSingle(Acquisition::Seismogram<ValueType>& seismo, lama::DenseVector<ValueType>& wavefieldSingle, lama::DenseVector<ValueType>& temp, IndexType t);

                /* wavefields */
                Wavefields::Wavefields<ValueType>& wavefield;
                
                /* source */
                Acquisition::SeismogramHandler<ValueType> const& sources;
                
                /* receiver */
                Acquisition::SeismogramHandler<ValueType>& receivers;
                
                /* Temporary memory */
                lama::DenseVector<ValueType> applySource_samplesVX;
                lama::DenseVector<ValueType> applySource_samplesVY;
                lama::DenseVector<ValueType> applySource_samplesVZ;
                
                lama::DenseVector<ValueType> gatherSeismogram_samplesVX;
                lama::DenseVector<ValueType> gatherSeismogram_samplesVY;
                lama::DenseVector<ValueType> gatherSeismogram_samplesVZ;

            };
            
        }
    }
}

template<typename ValueType>
KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImpl(Acquisition::Sources<ValueType> const& sourcesIN,Acquisition::Receivers<ValueType>& receiversIN, Wavefields::Wavefields<ValueType>& wavefieldIN)
:wavefield(wavefieldIN),sources(sourcesIN.getSeismogramHandler()),receivers(receiversIN.getSeismogramHandler())
{
    
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogram(IndexType t)
{
    if(receivers.getNumTracesGlobal(Acquisition::SeismogramType::P)>0){
        gatherSeismogramPressure(receivers.getSeismogram(Acquisition::SeismogramType::P),wavefield,t);
    }
    if(receivers.getNumTracesGlobal(Acquisition::SeismogramType::VX)>0){
        gatherSeismogramVX(receivers.getSeismogram(Acquisition::SeismogramType::VX),wavefield,t);
    }
    if(receivers.getNumTracesGlobal(Acquisition::SeismogramType::VY)>0){
        gatherSeismogramVY(receivers.getSeismogram(Acquisition::SeismogramType::VY),wavefield,t);
    }
    if(receivers.getNumTracesGlobal(Acquisition::SeismogramType::VZ)>0){
        gatherSeismogramVZ(receivers.getSeismogram(Acquisition::SeismogramType::VZ),wavefield,t);
    }
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogramVX(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    gatherSeismogramSingle(seismo,wavefieldIN.getVX(),gatherSeismogram_samplesVX,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogramVY(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    gatherSeismogramSingle(seismo,wavefieldIN.getVY(),gatherSeismogram_samplesVY,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogramVZ(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    gatherSeismogramSingle(seismo,wavefieldIN.getVZ(),gatherSeismogram_samplesVZ,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogramSingle(Acquisition::Seismogram<ValueType>& seismo, lama::DenseVector<ValueType>& wavefieldSingle, lama::DenseVector<ValueType>& temp, IndexType t)
{
    const lama::DenseVector<IndexType>& coordinates=seismo.getCoordinates();
    lama::DenseMatrix<ValueType>& seismogramData=seismo.getData();
    
    temp.gather(wavefieldSingle,coordinates,utilskernel::binary::BinaryOp::COPY);
    seismogramData.setColumn(temp,t,utilskernel::binary::BinaryOp::COPY);
}


template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySource(IndexType t)
{
    if(sources.getNumTracesGlobal(Acquisition::SeismogramType::P)>0){
        applySourcePressure(sources.getSeismogram(Acquisition::SeismogramType::P),wavefield,t);
    }
    if(sources.getNumTracesGlobal(Acquisition::SeismogramType::VX)>0){
        applySourceVX(sources.getSeismogram(Acquisition::SeismogramType::VX),wavefield,t);
    }
    if(sources.getNumTracesGlobal(Acquisition::SeismogramType::VY)>0){
        applySourceVY(sources.getSeismogram(Acquisition::SeismogramType::VY),wavefield,t);
    }
    if(sources.getNumTracesGlobal(Acquisition::SeismogramType::VZ)>0){
        applySourceVZ(sources.getSeismogram(Acquisition::SeismogramType::VZ),wavefield,t);
    }
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVX(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    applySourceSingle(seismo,wavefieldIN.getVX(),applySource_samplesVX,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVY(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    applySourceSingle(seismo,wavefieldIN.getVY(),applySource_samplesVY,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVZ(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t)
{
    applySourceSingle(seismo,wavefieldIN.getVZ(),applySource_samplesVZ,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceSingle(Acquisition::Seismogram<ValueType> const& seismo, lama::DenseVector<ValueType>& wavefieldSingle, lama::DenseVector<ValueType>& temp, IndexType t)
{
    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType>& sourcesSignals=seismo.getData();
    const lama::DenseVector<IndexType>& coordinates=seismo.getCoordinates();
    
    sourcesSignals.getColumn(temp,t);
    wavefieldSingle.scatter(coordinates,temp,utilskernel::binary::BinaryOp::ADD);
}
