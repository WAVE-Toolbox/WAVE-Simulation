#include "SourceReceiverImpl.hpp"
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
    setContextPtrToTemporary(wavefield.getContextPtr());
}

/*! \brief Gather seismogram.
 *
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogram(IndexType t)
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
    gatherSeismogramSingle(seismo, wavefieldIN.getRefEX(), gatherSeismogram_samplesEX, t);
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
    gatherSeismogramSingle(seismo, wavefieldIN.getRefEY(), gatherSeismogram_samplesEY, t);
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
    gatherSeismogramSingle(seismo, wavefieldIN.getRefEZ(), gatherSeismogram_samplesEZ, t);
}

/*! \brief Gather Hz seismogram.
 *
 \param seismo Seismogram
 \param wavefieldIN incoming wavefield
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramHZ(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    gatherSeismogramSingle(seismo, wavefieldIN.getRefHZ(), gatherSeismogram_samplesHZ, t);
}

/*! \brief Gather single seismogram.
 *
 \param seismo Seismogram
 \param wavefieldSingle Wavefield-single
  \param temp temporary Value
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::gatherSeismogramSingle(Acquisition::Seismogram<ValueType> &seismo, lama::DenseVector<ValueType> &wavefieldSingle, lama::DenseVector<ValueType> &temp, scai::IndexType t)
{
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();
    lama::DenseMatrix<ValueType> &seismogramData = seismo.getData();

    temp.gatherInto(wavefieldSingle, coordinates, common::BinaryOp::COPY);
    seismogramData.setColumn(temp, t, common::BinaryOp::COPY);
}

/*! \brief Applying source.
 *
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySource(IndexType t)
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
    applySourceSingle(seismo, wavefieldIN.getRefEX(), applySource_samplesEX, t);
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
    applySourceSingle(seismo, wavefieldIN.getRefEY(), applySource_samplesEY, t);
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
    applySourceSingle(seismo, wavefieldIN.getRefEZ(), applySource_samplesEZ, t);
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
    applySourceSingle(seismo, wavefieldIN.getRefHZ(), applySource_samplesHZ, t);
}

/*! \brief Applying single source.
 *
 \param seismo Seismogram
 \param wavefieldSingle Wavefields
 \param temp temporary Value
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<ValueType>::applySourceSingle(Acquisition::Seismogram<ValueType> const &seismo, lama::DenseVector<ValueType> &wavefieldSingle, lama::DenseVector<ValueType> &temp, scai::IndexType t)
{
    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignals = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();

    sourcesSignals.getColumn(temp, t);
    wavefieldSingle.scatter(coordinates, true, temp, common::BinaryOp::ADD);
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

template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<float>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImplEM<double>;
