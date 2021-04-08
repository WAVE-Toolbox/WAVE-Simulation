#include "FDTD2Delastic.hpp"
using namespace scai;

/*! \brief Gather the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getRefSyy();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();

    gatherSeismogram_samplesPressure.gatherInto(Sxx, coordinates, common::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gatherInto(Syy, coordinates, common::BinaryOp::ADD);
    gatherSeismogram_samplesPressure*=0.5;
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
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefield.getRefSyy();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.get1DCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    Sxx.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
    Syy.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD2Delastic<float>;
