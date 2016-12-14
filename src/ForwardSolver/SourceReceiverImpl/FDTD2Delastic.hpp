
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
        
        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl {
            
            template<typename ValueType>
            class FDTD2Delastic : public SourceReceiverImpl<ValueType>
            {
            public:
                
                FDTD2Delastic()=delete;
                ~FDTD2Delastic(){};
                
                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;
                
                inline void applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t) override;
                inline void gatherSeismogramPressure(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t) override;

            private:
                
                /* Temporary memory */
                lama::DenseVector<ValueType> applySource_samplesPressure;
                lama::DenseVector<ValueType> gatherSeismogram_samplesPressure;
            };
            
        }
    }
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
    lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
    
    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType>& seismogramDataPressure=seismo.getData();
    const lama::DenseVector<IndexType>& coordinates=seismo.getCoordinates();

    gatherSeismogram_samplesPressure.gather(Sxx,coordinates,utilskernel::binary::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gather(Syy,coordinates,utilskernel::binary::BinaryOp::ADD);
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure,t,utilskernel::binary::BinaryOp::COPY);
}

template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
    
        /* Get reference to sourcesignal storing seismogram */
        const lama::DenseMatrix<ValueType>& sourcesSignalsPressure=seismo.getData();
        const lama::DenseVector<IndexType>& coordinatesPressure=seismo.getCoordinates();
        
        sourcesSignalsPressure.getColumn(applySource_samplesPressure,t);
        Sxx.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
        Syy.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
    
}

