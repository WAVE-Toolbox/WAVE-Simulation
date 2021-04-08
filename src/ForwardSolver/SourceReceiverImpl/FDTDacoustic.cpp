#include "FDTDacoustic.hpp"
using namespace scai;

/*! \brief Gather the seismogram pressure.
 *
 *
 \param seismo Seismogram
 \param wavefield Wavefields
 \param t Time-step
 */
template <typename ValueType>
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<ValueType>::gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &p = wavefield.getRefP();

    /* Gather seismogram for the pressure traces */
    lama::DenseMatrix<ValueType> &seismogramDataPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinates = seismo.get1DCoordinates();

    gatherSeismogram_samplesPressure.gatherInto(p, coordinates, common::BinaryOp::COPY);
    gatherSeismogram_samplesPressure*=1;
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
void KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<ValueType>::applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, scai::IndexType t)
{
    /* Get reference to wavefields */
    lama::DenseVector<ValueType> &p = wavefield.getRefP();

    /* Get reference to sourcesignal storing seismogram */
    const lama::DenseMatrix<ValueType> &sourcesSignalsPressure = seismo.getData();
    const lama::DenseVector<IndexType> &coordinatesPressure = seismo.get1DCoordinates();

    sourcesSignalsPressure.getColumn(applySource_samplesPressure, t);
    p.scatter(coordinatesPressure, true, applySource_samplesPressure, common::BinaryOp::ADD);
}

template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::FDTDacoustic<float>;
