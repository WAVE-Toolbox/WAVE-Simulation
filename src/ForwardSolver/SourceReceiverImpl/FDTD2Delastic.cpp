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
    lama::DenseVector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getRefSyy();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.getCoordinates();

    gatherSeismogram_samplesPressure.gather(Sxx, coordinates, common::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gather(Syy, coordinates, common::BinaryOp::ADD);
    gatherSeismogram_samplesPressure*=-0.5;
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure, t, common::BinaryOp::COPY);
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
    lama::DenseVector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getRefSyy();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.getCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    Sxx.scatter(coordinatesPressure, applySource_samplesPressure, common::BinaryOp::ADD);
    Syy.scatter(coordinatesPressure, applySource_samplesPressure, common::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<float>;
