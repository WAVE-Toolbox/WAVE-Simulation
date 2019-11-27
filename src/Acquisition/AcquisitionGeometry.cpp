#include "AcquisitionGeometry.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initOptionalAcquisitionParameter(IndexType /*numParameter*/, IndexType /*numTracesGlobal*/, scai::lama::DenseMatrix<ValueType> /*acquisition*/, scai::dmemo::DistributionPtr /*dist_wavefield_traces*/, hmemo::ContextPtr /*ctx*/)
{
}

/*! \brief Initialization of the seismogram handler
 *
 *
 \param NT Number of timesteps
 \param ctx Context
 \param dist_wavefield Distribution of global grid
 */
template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initSeismogramHandler(IndexType const NT, scai::hmemo::ContextPtr const ctx, scai::dmemo::DistributionPtr const dist_wavefield)
{

    SCAI_ASSERT_DEBUG(seismogramTypes.size() == coordinates1D.size(), "Size mismatch");
    SCAI_ASSERT_DEBUG(numTracesGlobal == seismogramTypes.size(), "Size mismatch");

    IndexType count[NUM_ELEMENTS_SEISMOGRAMTYPE] = {0, 0, 0, 0};
    lama::DenseVector<IndexType> coord[NUM_ELEMENTS_SEISMOGRAMTYPE];
    dmemo::DistributionPtr dist[NUM_ELEMENTS_SEISMOGRAMTYPE];

    /* Count elements for each source type */
    IndexType tempIndexType;
    for (IndexType i = 0; i < numTracesGlobal; ++i) {
        tempIndexType = seismogramTypes[i] - 1;
        SCAI_ASSERT_VALID_INDEX_DEBUG(tempIndexType, NUM_ELEMENTS_SEISMOGRAMTYPE, "Unknown Type in trace Nr. " << i + 1);

        ++count[tempIndexType];
    }

    /* Allocate lama vectors */
    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; ++i) {
        coord[i].allocate(count[i]);
        count[i] = 0;
    }

    /* Sort coordinates */
    for (IndexType i = 0; i < numTracesGlobal; ++i) {

        tempIndexType = seismogramTypes[i] - 1;

        coord[tempIndexType].setValue(count[tempIndexType], this->get1DCoordinates().getValue(i));
        ++count[tempIndexType];
    }

    SCAI_ASSERT_EQ_DEBUG(static_cast<SeismogramType>(0), SeismogramType::P, "Cast went wrong");
    SCAI_ASSERT_EQ_DEBUG(static_cast<SeismogramType>(1), SeismogramType::VX, "Cast went wrong");
    SCAI_ASSERT_EQ_DEBUG(static_cast<SeismogramType>(2), SeismogramType::VY, "Cast went wrong");
    SCAI_ASSERT_EQ_DEBUG(static_cast<SeismogramType>(3), SeismogramType::VZ, "Cast went wrong");

    /* Calculate distribution, redistribute coordinates and set coordinates to seismogramHandler */
    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; ++i) {
        if (coord[i].size() > 0) {
            dist[i] = calcDistribution(coord[i], dist_wavefield);
            seismograms.getSeismogram(static_cast<SeismogramType>(i)).allocate(ctx, dist[i], NT);
            coord[i].redistribute(dist[i]);
            seismograms.getSeismogram(static_cast<SeismogramType>(i)).setCoordinates(coord[i]);
        }
        count[i] = 0;
    }

    seismograms.setContextPtr(ctx);
}

/*! \brief get number of parameters
 *
 * returns the number of parameters of the aquisition matrix eg. X,Y,Z,Seismogram Type
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getNumParameter() const
{
    return numParameter;
}


/*! \brief Get reference to the #SeismogramHandler
 *
 *
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType> const &KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getSeismogramHandler() const
{
    return (seismograms);
}

template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType> &KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getSeismogramHandler()
{
    return (seismograms);
}

/*! \brief Get reference to traces type as #SeismogramType
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getSeismogramTypes() const
{
    SCAI_ASSERT_DEBUG(numTracesGlobal == seismogramTypes.size(), "Size mismatch");
    return (seismogramTypes);
}

/*! \brief Get reference to the coordinates
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::AcquisitionGeometry<ValueType>::get1DCoordinates() const
{
    SCAI_ASSERT_DEBUG(numTracesGlobal == coordinates1D.size(), "Size mismatch");
    return (coordinates1D);
}

/*! \brief Get number of global traces
 *
 \return Number of global traces
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getNumTracesGlobal() const
{
    return (numTracesGlobal);
}

/*! \brief Get number of local traces
 *
 \return Number of local traces on this process
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getNumTracesLocal() const
{
    return (numTracesLocal);
}

/*! \brief Write acquisition to file
 *
 \param filename Filename to write acquisition
 */
template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::writeAcquisitionToFile(std::string const & /*filename*/) const
{
}

template class KITGPI::Acquisition::AcquisitionGeometry<double>;
template class KITGPI::Acquisition::AcquisitionGeometry<float>;
