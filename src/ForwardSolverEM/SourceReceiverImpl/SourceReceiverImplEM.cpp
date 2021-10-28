#include "SourceReceiverImplEM.hpp"
using namespace scai;

/*! \brief default method for Source-Receiver Implementation.
 *
 \param sourceConfig Sources
 \param receiverConfig Recievers
 \param wavefieldIN Wavefield
 */
template <typename ValueType>
KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::SourceReceiverImplEM(Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN)
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
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogram(scai::IndexType t)
{
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EZ) > 0) {
        gatherSeismogramEZ(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramTypeEM::EZ), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EX) > 0) {
        gatherSeismogramEX(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramTypeEM::EX), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EY) > 0) {
        gatherSeismogramEY(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramTypeEM::EY), wavefield, t);
    }
    if (seismogramHandlerRec.getNumTracesGlobal(Acquisition::SeismogramTypeEM::HZ) > 0) {
        gatherSeismogramHZ(seismogramHandlerRec.getSeismogram(Acquisition::SeismogramTypeEM::HZ), wavefield, t);
    }
}

/*! \brief Gather Ex seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramEX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefEX(), gatherSeismogram_samplesEX, t);
}

/*! \brief Gather Ey seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramEY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefEY(), gatherSeismogram_samplesEY, t);
}

/*! \brief Gather Ez seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramEZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefEZ(), gatherSeismogram_samplesEZ, t);
}

/*! \brief Gather Hz seismogram.
 *
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->gatherSeismogramSingle(seismo, wavefieldIN.getRefHZ(), gatherSeismogram_samplesHZ, t);
}

/*! \brief Applying source.
 *
 \param seismogramHandlerSrc Seismogram handler source
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySource(scai::IndexType t)
{
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EZ) > 0) {
        applySourceEZ(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramTypeEM::EZ), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EX) > 0) {
        applySourceEX(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramTypeEM::EX), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramTypeEM::EY) > 0) {
        applySourceEY(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramTypeEM::EY), wavefield, t);
    }
    if (seismogramHandlerSrc.getNumTracesGlobal(Acquisition::SeismogramTypeEM::HZ) > 0) {
        applySourceHZ(seismogramHandlerSrc.getSeismogram(Acquisition::SeismogramTypeEM::HZ), wavefield, t);
    }
}

/*! \brief Applying Ex source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceEX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefEX(), applySource_samplesEX, t);
}

/*! \brief Applying Ey source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceEY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefEY(), applySource_samplesEY, t);
}

/*! \brief Applying Ez source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceEZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefEZ(), applySource_samplesEZ, t);
}

/*! \brief Applying Hz source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceHZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    this->applySourceSingle(seismo, wavefieldIN.getRefHZ(), applySource_samplesHZ, t);
}

/*! \brief Setting context pointer to temporary source and seismogram.
 *
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::setContextPtrToTemporary(hmemo::ContextPtr ctx)
{
    applySource_samplesEX.setContextPtr(ctx);
    applySource_samplesEY.setContextPtr(ctx);
    applySource_samplesEZ.setContextPtr(ctx);
    applySource_samplesHZ.setContextPtr(ctx);
    gatherSeismogram_samplesEX.setContextPtr(ctx);
    gatherSeismogram_samplesEY.setContextPtr(ctx);
    gatherSeismogram_samplesEZ.setContextPtr(ctx);
    gatherSeismogram_samplesHZ.setContextPtr(ctx);
}



/* Seismic */
/*! \brief Gather the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefieldIN Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramPressure in an EM modelling")
}

/*! \brief Gather vx seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramVX(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramVX in an EM modelling")
}

/*! \brief Gather vy seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramVY(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramVY in an EM modelling")
}

/*! \brief Gather vz seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no gatherSeismogramVZ in an EM modelling")
}

/*! \brief Applying pressure from source.
 *
 *
 \param seismo Seismogram
 \param wavefieldIN Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourcePressure in an EM modelling")
}

/*! \brief Applying vx source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceVX(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceVX in an EM modelling")
}

/*! \brief Applying vy source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceVY(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceVY in an EM modelling")
}

/*! \brief Applying vz source.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming Wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceVZ(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    COMMON_THROWEXCEPTION("There is no applySourceVZ in an EM modelling")
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<float>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<double>;
