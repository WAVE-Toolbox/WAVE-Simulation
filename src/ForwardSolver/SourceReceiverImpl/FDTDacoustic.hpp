
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
        
        namespace SourceReceiverImpl {
            
            //! \brief FDTDacoustic class
            template<typename ValueType>
            class FDTDacoustic : public SourceReceiverImpl<ValueType>
            {
            public:
                
                //! Default constructor
                FDTDacoustic()=delete;
                //! Default destructor
                ~FDTDacoustic(){};
                
                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;
                
                inline void applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t) override;
                inline void gatherSeismogramPressure(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t) override;

            private:
                
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;

            };
            
        }
    }
}

/*! \brief Gether the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType>& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType>& p=wavefield.getP();
    
    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType>& seismogramDataPressure=seismo.getData();
    const lama::DenseVector<IndexType>& coordinates=seismo.getCoordinates();

    gatherSeismogram_samplesPressure.gather(p,coordinates,utilskernel::binary::BinaryOp::COPY);
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure,t,utilskernel::binary::BinaryOp::COPY);
}

/*! \brief Applying pressure from source.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template<typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const& seismo, Wavefields::Wavefields<ValueType>& wavefield, IndexType t)
{
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& p=wavefield.getP();
    
        /* Get reference to sourcesignal storing seismogram */
        const lama::DenseMatrix<ValueType>& sourcesSignalsPressure=seismo.getData();
        const lama::DenseVector<IndexType>& coordinatesPressure=seismo.getCoordinates();
        
        sourcesSignalsPressure.getColumn(applySource_samplesPressure,t);
        p.scatter(coordinatesPressure,applySource_samplesPressure,utilskernel::binary::BinaryOp::ADD);
    
}

