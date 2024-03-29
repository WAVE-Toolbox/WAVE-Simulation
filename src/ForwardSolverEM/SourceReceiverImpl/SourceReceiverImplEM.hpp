#pragma once

#include <scai/lama.hpp>

#include "../../ForwardSolver/SourceReceiverImpl/SourceReceiverImpl.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief SourceReceiverImpl class
            template <typename ValueType>
            class SourceReceiverImplEM : public SourceReceiverImpl<ValueType>
            {
              public:
                //! Default constructor
                SourceReceiverImplEM() = delete;
                explicit SourceReceiverImplEM(Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN);
                //! Default destructor
                ~SourceReceiverImplEM(){};
                           
                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;
                
                void applySource(scai::IndexType t) override;
                void gatherSeismogram(scai::IndexType t) override;

              protected:           
                /* Common */  
                void setContextPtrToTemporary(scai::hmemo::ContextPtr ctx) override; 
                /* wavefields */
                Wavefields::Wavefields<ValueType> &wavefield; //!< Wavefields

                /* source */
                Acquisition::SeismogramHandler<ValueType> const &seismogramHandlerSrc; //!< Seismogram Handler of Sources

                /* receiver */
                Acquisition::SeismogramHandler<ValueType> &seismogramHandlerRec; //!< Seismogram Handler of Receivers

                /* Seismic */
                //! \brief Allying method for pressure source
                virtual void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceVX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceVY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceVZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramVX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramVY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesVX; //!< Variable to apply velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::applySource_samplesVY; //!< Variable to apply velocity Source-samples in y-direction
                using SourceReceiverImpl<ValueType>::applySource_samplesVZ; //!< Variable to apply velocity Source-samples in z-direction

                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesVX; //!< Variable to gather velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesVY; //!< Variable to gather velocity Source-samples in y-direction
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesVZ; //!< Variable to gather velocity Source-samples in z-direction

                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;      //!< Variable to apply pressure Source-samples
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure; //!< Variable to gather pressure Source-samples

                /* EM */
                //! \brief Allying method for pressure source
                virtual void applySourceEZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceEX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceEY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void applySourceHZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramEX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramEY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramEZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                virtual void gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesEZ; //!< Variable to apply velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::applySource_samplesEX; //!< Variable to apply velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::applySource_samplesEY; //!< Variable to apply velocity Source-samples in y-direction
                using SourceReceiverImpl<ValueType>::applySource_samplesHZ; //!< Variable to apply velocity Source-samples in z-direction

                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesEZ; //!< Variable to gather velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesEX; //!< Variable to gather velocity Source-samples in x-direction
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesEY; //!< Variable to gather velocity Source-samples in y-direction
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesHZ; //!< Variable to gather velocity Source-samples in z-direction
                
            };
        }
    }
}
