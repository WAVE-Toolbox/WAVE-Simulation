#pragma once

#include <scai/lama.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/AcquisitionGeometry.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief SourceReceiverImpl class
            template <typename ValueType>
            class SourceReceiverImpl
            {
              public:
                //! Default constructor
                SourceReceiverImpl() = delete;
                explicit SourceReceiverImpl(Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN);
                //! Default destructor
                ~SourceReceiverImpl(){};

                void applySource(scai::IndexType t);
                void gatherSeismogram(scai::IndexType t);

              protected:
                //! \brief Allying method for pressure source
                virtual void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceVX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void applySourceVY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void applySourceVZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                void applySourceSingle(Acquisition::Seismogram<ValueType> const &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramVX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void gatherSeismogramVY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                virtual void gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t);
                void gatherSeismogramSingle(Acquisition::Seismogram<ValueType> &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);

                /* wavefields */
                Wavefields::Wavefields<ValueType> &wavefield; //!< Wavefields

                /* source */
                Acquisition::SeismogramHandler<ValueType> const &seismogramHandlerSrc; //!< Seismogram Handler of Sources

                /* receiver */
                Acquisition::SeismogramHandler<ValueType> &seismogramHandlerRec; //!< Seismogram Handler of Receivers

                /* Temporary memory */
                scai::lama::DenseVector<ValueType> applySource_samplesVX; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesVY; //!< Variable to apply velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> applySource_samplesVZ; //!< Variable to apply velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVX; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVY; //!< Variable to gather velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVZ; //!< Variable to gather velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> applySource_samplesPressure;      //!< Variable to apply pressure Source-samples
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesPressure; //!< Variable to gather pressure Source-samples

              private:
                void setContextPtrToTemporary(scai::hmemo::ContextPtr ctx);
            };
        }
    }
}
