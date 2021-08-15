#include "SourceReceiverImpl.hpp"
using namespace scai;

/*! \brief Gather single seismogram.
 *
 \param seismo Seismogram
 \param wavefieldSingle Wavefield-single
  \param temp temporary Value
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::gatherSeismogramSingle(Acquisition::Seismogram<ValueType> &seismo, lama::DenseVector<ValueType> &wavefieldSingle, lama::DenseVector<ValueType> &temp, scai::IndexType t)
{
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();
    lama::DenseMatrix<ValueType> &seismogramData = seismo.getData();

    temp.gatherInto(wavefieldSingle, coordinates, common::BinaryOp::COPY);
    seismogramData.setColumn(temp, t, common::BinaryOp::COPY);
}

/*! \brief Applying single source.
 *
 \param seismo Seismogram
 \param wavefieldSingle Wavefields
 \param temp temporary Value
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::applySourceSingle(Acquisition::Seismogram<ValueType> const &seismo, lama::DenseVector<ValueType> &wavefieldSingle, lama::DenseVector<ValueType> &temp, scai::IndexType t)
{
    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignals = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();

    sourcesSignals.getColumn(temp, t);
    wavefieldSingle.scatter(coordinates, true, temp, common::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<float>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<double>;
