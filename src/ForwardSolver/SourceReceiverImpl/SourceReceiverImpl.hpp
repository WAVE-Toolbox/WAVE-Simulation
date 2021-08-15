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
                //! \brief Declare SourceReceiverImpl pointer
                typedef std::shared_ptr<SourceReceiverImpl<ValueType>> SourceReceiverImplPtr;
                
                //! Default constructor
                SourceReceiverImpl(){};
                //! Default destructor
                ~SourceReceiverImpl(){};
                                
                virtual void init(Acquisition::SeismogramHandler<ValueType> const &seismogramHandlerSrc, Acquisition::SeismogramHandler<ValueType> &seismogramHandlerRec, Wavefields::Wavefields<ValueType> &wavefieldIN) = 0;
                virtual void applySource(Acquisition::SeismogramHandler<ValueType> const &seismogramHandlerSrc, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogram(Acquisition::SeismogramHandler<ValueType> &seismogramHandlerRec, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;

              protected:                
                /* Common */
                void applySourceSingle(Acquisition::Seismogram<ValueType> const &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);
                void gatherSeismogramSingle(Acquisition::Seismogram<ValueType> &seismo, scai::lama::DenseVector<ValueType> &wavefieldSingle, scai::lama::DenseVector<ValueType> &temp, scai::IndexType t);
                virtual void setContextPtrToTemporary(scai::hmemo::ContextPtr ctx) = 0;
                
                /* Seismic */
                //! \brief Allying method for pressure source
                virtual void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceVX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceVY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceVZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramVX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramVY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                
                /* Temporary memory */
                scai::lama::DenseVector<ValueType> applySource_samplesVX; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesVY; //!< Variable to apply velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> applySource_samplesVZ; //!< Variable to apply velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVX; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVY; //!< Variable to gather velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesVZ; //!< Variable to gather velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> applySource_samplesPressure;      //!< Variable to apply pressure Source-samples
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesPressure; //!< Variable to gather pressure Source-samples

                /* EM */
                //! \brief Allying method for pressure source
                virtual void applySourceEZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceEX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceEY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void applySourceHZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;

                //! \brief Gather method for pressure seismogram
                virtual void gatherSeismogramEX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramEY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramEZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                virtual void gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) = 0;
                
                /* Temporary memory */
                scai::lama::DenseVector<ValueType> applySource_samplesEZ; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesEX; //!< Variable to apply velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> applySource_samplesEY; //!< Variable to apply velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> applySource_samplesHZ; //!< Variable to apply velocity Source-samples in z-direction

                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEZ; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEX; //!< Variable to gather velocity Source-samples in x-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesEY; //!< Variable to gather velocity Source-samples in y-direction
                scai::lama::DenseVector<ValueType> gatherSeismogram_samplesHZ; //!< Variable to gather velocity Source-samples in z-direction                
            };
        }
    }
}
