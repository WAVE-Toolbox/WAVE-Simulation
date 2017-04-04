#include "FDTD2Delastic.hpp"
using namespace scai;

/*! \brief Gether the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefield.getSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getSyy();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.getCoordinates();

    gatherSeismogram_samplesPressure.gather(Sxx, coordinates, scai::utilskernel::binary::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gather(Syy, coordinates, scai::utilskernel::binary::BinaryOp::ADD);
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure, t, scai::utilskernel::binary::BinaryOp::COPY);
}

/*! \brief Applying pressure from source.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefield.getSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getSyy();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.getCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    Sxx.scatter(coordinatesPressure, applySource_samplesPressure, scai::utilskernel::binary::BinaryOp::ADD);
    Syy.scatter(coordinatesPressure, applySource_samplesPressure, scai::utilskernel::binary::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<float>;
