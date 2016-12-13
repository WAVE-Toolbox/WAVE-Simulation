
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
                explicit SourceReceiverImpl(Acquisition::Sources<ValueType> const& sourcesIN, Wavefields::Wavefields<ValueType>& wavefieldIN);
                ~SourceReceiverImpl(){};
                
                void applySource(IndexType t);
                
            protected:
                
                virtual void applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) = 0;
                
                virtual void applySourceVX(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;
                virtual void applySourceVY(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;
                virtual void applySourceVZ(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefieldIN, IndexType t) ;

                inline void applySourceSingle(Acquisition::Seismogram<ValueType> const& seismo, lama::DenseVector<ValueType>& wavefieldSingle, lama::DenseVector<ValueType>& temp, IndexType t);
                
                Wavefields::Wavefields<ValueType>& wavefield;
                
                /* source */
                Acquisition::SeismogramHandler<ValueType> const& sources;
                
                /* receiver */
                //Acquisition::SeismogramHandler<ValueType> const& seismograms;
                
                
                /* Temporary memory */
                lama::DenseVector<ValueType> applySource_samplesVX;
                lama::DenseVector<ValueType> applySource_samplesVY;
                lama::DenseVector<ValueType> applySource_samplesVZ;
            };
            
        }
    }
}

template<typename ValueType>
KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImpl(Acquisition::Sources<ValueType> const& sourcesIN, Wavefields::Wavefields<ValueType>& wavefieldIN)
:wavefield(wavefieldIN),sources(sourcesIN.getSeismogramHandler())
{
    
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
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVX(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
    applySourceSingle(seismo,wavefield.getVX(),applySource_samplesVX,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVY(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
    applySourceSingle(seismo,wavefield.getVY(),applySource_samplesVY,t);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceVZ(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
    applySourceSingle(seismo,wavefield.getVZ(),applySource_samplesVZ,t);
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
