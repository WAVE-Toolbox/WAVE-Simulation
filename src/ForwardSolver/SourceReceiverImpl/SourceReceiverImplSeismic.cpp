#include "SourceReceiverImplSeismic.hpp"
using namespace scai;

/*! \brief default method for Source-Receiver Implementation.
 *
 \param sourceConfig Sources
 \param receiverConfig Recievers
 \param wavefieldIN Wavefield
 */
template <typename ValueType>
KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::SourceReceiverImplSeismic(Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN)
 : wavefield(wavefieldIN), seismogramHandlerSrc(sourceConfig.getSeismogramHandler()), seismogramHandlerRec(receiverConfig.getSeismogramHandler())
{
    // Set source coordinate to receiver seismogram handler
    if (seismogramHandlerSrc.getNumTracesTotal() == 1) {
        /* If only one source is injected in this simulation, the coordinate of this source is
         set to the receiver seismogram handler */
        IndexType temp = sourceConfig.get1DCoordinates().getValue(0);
        seismogramHandlerRec.setSourceCoordinate(temp);
    } else {
        /* If more than one source is injected at the same time, the source coordinate is
        set to zero, due to the choice of one coordinate would be subjective */
        seismogramHandlerRec.setSourceCoordinate(0);
    }
    // Set Context of temporary variables to context of wavefields
    setContextPtrToTemporary(wavefieldIN.getContextPtr());
}

/*! \brief Gather seismogram.
 *
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogram(scai::IndexType t)
{
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramType::P) > 0) {
        gatherSeismogramPressure(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramType::P), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramType::VX) > 0) {
        gatherSeismogramVX(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramType::VX), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramType::VY) > 0) {
        gatherSeismogramVY(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramType::VY), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramType::VZ) > 0) {
        gatherSeismogramVZ(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramType::VZ), wavefield, t);
    }
}

/*! \brief Gather vx seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramVX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefVX(), gatherSeismogram_samplesVX, t);
}

/*! \brief Gather vy seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramVY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefVY(), gatherSeismogram_samplesVY, t);
}

/*! \brief Gather vz seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefVZ(), gatherSeismogram_samplesVZ, t);
}

/*! \brief Applying source.
 *
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySource(scai::IndexType t)
{
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramType::P) > 0) {
        applySourcePressure(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramType::P), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramType::VX) > 0) {
        applySourceVX(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramType::VX), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramType::VY) > 0) {
        applySourceVY(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramType::VY), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramType::VZ) > 0) {
        applySourceVZ(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramType::VZ), wavefield, t);
    }
}

/*! \brief Applying vx source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceVX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefVX(), applySource_samplesVX, t);
}

/*! \brief Applying vy source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceVY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefVY(), applySource_samplesVY, t);
}

/*! \brief Applying vz source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceVZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefVZ(), applySource_samplesVZ, t);
}

/*! \brief Setting context pointer to temporary source and seismogram.
 *
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::setContextPtrToTemporary(hmemo::ContextPtr ctx)
{
    applySource_samplesVX.setContextPtr(ctx);
    applySource_samplesVY.setContextPtr(ctx);
    applySource_samplesVZ.setContextPtr(ctx);
    gatherSeismogram_samplesVX.setContextPtr(ctx);
    gatherSeismogram_samplesVY.setContextPtr(ctx);
    gatherSeismogram_samplesVZ.setContextPtr(ctx);
    applySource_samplesPressure.setContextPtr(ctx);
    gatherSeismogram_samplesPressure.setContextPtr(ctx);
}



/* EM */
/*! \brief Gather Ex seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramEX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramEX in an Seismic modelling")
}

/*! \brief Gather Ey seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramEY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramEY in an Seismic modelling")
}

/*! \brief Gather Ez seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramEZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramEZ in an Seismic modelling")
}

/*! \brief Gather Hz seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramHZ in an Seismic modelling")
}

/*! \brief Applying Ex source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceEX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceEX in an Seismic modelling")
}

/*! \brief Applying Ey source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceEY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceEY in an Seismic modelling")
}

/*! \brief Applying Ez source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceEZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceEZ in an Seismic modelling")
}

/*! \brief Applying Hz source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<ValueType>::applySourceHZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceHZ in an Seismic modelling")
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<float>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplSeismic<double>;
