#pragma once

#include <scai/lama.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/AcquisitionGeometry.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../WavefieldsEM/Wavefields.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief SourceReceiverImpl class
            template <typename ValueType>
            class SourceReceiverImplEM
            {
              public:
                //! \brief Declare SourceReceiverImpl pointer
                typedef std::shared_ptr<SourceReceiverImplEM<ValueType>> SourceReceiverImplPtr;
                
                //! Default constructor
                SourceReceiverImplEM() = delete;
                explicit SourceReceiverImplEM(Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::WavefieldsEM<ValueType> &wavefieldIN);
                //! Default destructor
                ~SourceReceiverImplEM(){};

                void applySource(scai::IndexType t);
                void gatherSeismogram(scai::IndexType t);

              protected:
                //! \brief Allying method for pressure source
                virtual void applySourceEZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void applySourceEX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void applySourceEY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void applySourceHZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                void applySourceSingle(Acquisition::Seismogram<ValueType> const &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramEX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void gatherSeismogramEY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void gatherSeismogramEZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::WavefieldsEM<ValueType> &wavefieldIN, scai::IndexType t);
                void gatherSeismogramSingle(Acquisition::Seismogram<ValueType> &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);

                /* wavefields */
                Wavefields::WavefieldsEM<ValueType> &wavefield; //!< Wavefields

                /* source */
                Acquisition::SeismogramHandler<ValueType> const &seismogramHandlerSrc; //!< Seismogram Handler of Sources

                /* receiver */
                Acquisition::SeismogramHandler<ValueType> &seismogramHandlerRec; //!< Seismogram Handler of Receivers

                /* Temporary memory */
                scai::lama::DenseVector<ValueType> applySource_samplesEZ; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesEX; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesEY; //!< Variable to apply velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> applySource_samplesHZ; //!< Variable to apply velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEZ; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEX; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEY; //!< Variable to gather velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesHZ; //!< Variable to gather velocity Source-samples in z-direction

              private:
                void setContextPtrToTemporary(scai::hmemo::ContextPtr ctx);
            };
        }
    }
}
