#include "FDTD3Delastic.hpp"
using namespace scai;

/*! \brief Gather the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefieldIN Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD3Delastic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefieldIN.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefieldIN.getRefSyy();
    lama::DenseVector<ValueType> &Szz = wavefieldIN.getRefSzz();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();

    gatherSeismogram_samplesPressure.gatherInto(Sxx, coordinates, common::BinaryOp::COPY);
    gatherSeismogram_samplesPressure.gatherInto(Syy, coordinates, common::BinaryOp::ADD);
    gatherSeismogram_samplesPressure.gatherInto(Szz, coordinates, common::BinaryOp::ADD);
    gatherSeismogram_samplesPressure/=3;
    seismogramDataPressure.setColumn(gatherSeismogram_samplesPressure, t, common::BinaryOp::COPY);
}

/*! \brief Applying pressure from source.
 *
 *
 \param seismo Seismogram
 \param wavefieldIN Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTD3Delastic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &Sxx = wavefieldIN.getRefSxx();
    lama::DenseVector<ValueType> &Syy = wavefieldIN.getRefSyy();
    lama::DenseVector<ValueType> &Szz = wavefieldIN.getRefSzz();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.get1DCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    Sxx.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
    Syy.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
    Szz.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD3Delastic<float>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTD3Delastic<double>;
