#include "FDTD2Dsh.hpp"
using namespace scai;

/*! \brief Gether the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Dsh<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxz = wavefield.getSxz();
    lama::DenseVector<ValueType> &Syz = wavefield.getSyz();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.getCoordinates();

    gatherSeismogram_samplesPressure.gather(Sxz, coordinates, scai::common::binary::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gather(Syz, coordinates, scai::common::binary::BinaryOp::ADD);
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure, t, scai::common::binary::BinaryOp::COPY);
}

/*! \brief Applying pressure from source.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Dsh<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxz = wavefield.getSxz();
    lama::DenseVector<ValueType> &Syz = wavefield.getSyz();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.getCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    Sxz.scatter(coordinatesPressure, applySource_samplesPressure, scai::common::binary::BinaryOp::ADD);
    Syz.scatter(coordinatesPressure, applySource_samplesPressure, scai::common::binary::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Dsh<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Dsh<float>;
