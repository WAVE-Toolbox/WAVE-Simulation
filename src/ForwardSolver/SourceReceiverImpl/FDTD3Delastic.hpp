
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

#include "SourceReceiverImpl.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief Derivatives namespace
        namespace SourceReceiverImpl {
            
            template<typename ValueType>
            class FDTD3Delastic : public SourceReceiverImpl<ValueType>
            {
            public:
                
                FDTD3Delastic()=delete;
                ~FDTD3Delastic(){};
                
                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;
                
                void applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t) override;
                
            private:
                
                /* Temporary memory */
                lama::DenseVector<ValueType> applySource_samplesPressure;
            };
            
        }
    }
}


template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD3Delastic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
        lama::DenseVector<ValueType>& Szz=wavefield.getSzz();
        
        /* Get reference to sourcesignal storing seismogram */
        const lama::DenseMatrix<ValueType>& sourcesSignalsPressure=seismo.getData();
        const lama::DenseVector<IndexType>& coordinatesPressure=seismo.getCoordinates();
        
        sourcesSignalsPressure.getColumn(applySource_samplesPressure,t);
        Sxx.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
        Syy.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
        Szz.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
    
}

