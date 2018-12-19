#include "AcquisitionGeometry.hpp"

using namespace scai;

/*! \brief Determination of local indices based on given global indeces
 *
 * Calculate the number of indeces within the local processing unit as well as
 * the indeces of the local index.
 *
 \param coordinatesglobal DenseVector with global coordinates
 \param localIndices DenseVector with local coordinates
 \param dist Distribution of global grid
 */
template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::Global2Local(scai::lama::Vector<IndexType> const &coordinatesglobal, scai::hmemo::HArray<IndexType> &localIndices, scai::dmemo::DistributionPtr dist) const
{

    IndexType n_global = coordinatesglobal.size(); // Number of global entries

    IndexType coordinatetemp_int;

    IndexType i = 0;
    for (IndexType n = 0; n < n_global; n++) {

        coordinatetemp_int = coordinatesglobal.getValue(n);

        if (dist->isLocal(coordinatetemp_int)) {
            i++;
        }
    }

    /* Determine coordinates of local receivers in the global coordinate vector */
    localIndices.resize(i);
    hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
    i = 0;
    for (IndexType n = 0; n < n_global; n++) {

        coordinatetemp_int = coordinatesglobal.getValue(n);
        if (dist->isLocal(coordinatetemp_int)) {
            write_localIndices[i] = n;
            i++;
        }
    }
}

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initOptionalAcquisitionParameter(IndexType /*numParameter*/, IndexType /*numTracesGlobal*/, scai::lama::DenseMatrix<ValueType> /*acquisition*/, scai::dmemo::DistributionPtr /*dist_wavefield_traces*/, hmemo::ContextPtr /*ctx*/)
{
}

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initSeismogramHandler(IndexType const NT, scai::hmemo::ContextPtr const ctx, scai::dmemo::DistributionPtr const dist_wavefield)
{

    SCAI_ASSERT_DEBUG(seismogramTypes.size() == coordinates.size(), "Size mismatch");
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

        coord[tempIndexType].setValue(count[tempIndexType], this->getCoordinates().getValue(i));
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

template <typename ValueType>
IndexType KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getNumParameter() const
{
    return numParameter;
}

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::setAcquisition(scai::lama::DenseMatrix<ValueType> acquisition_temp, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx)
{

    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NX<=0");

    IndexType nrow_temp = acquisition_temp.getNumRows();
    IndexType ncolumn_temp = acquisition_temp.getNumColumns();

    /* Derive number of traces and number of read-in parameters */
    numTracesGlobal = nrow_temp;
    numParameter = ncolumn_temp;

    checkRequiredNumParameter(numParameter);

    /* Distribution: Master process only (comm->myRank()==0) */
    dmemo::DistributionPtr dist_master_numParameter(new dmemo::CyclicDistribution(numParameter, numParameter, dist_wavefield->getCommunicatorPtr()));
    dmemo::DistributionPtr dist_master_numTracesGlobal(new dmemo::CyclicDistribution(numTracesGlobal, numTracesGlobal, dist_wavefield->getCommunicatorPtr()));

    /* Distribution: Replicated on all processes */
    dmemo::DistributionPtr no_dist_numTracesGlobal(new scai::dmemo::NoDistribution(numTracesGlobal));
    dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(numParameter));

    /* Allocate acquisition matrix on master */
    acquisition.allocate(dist_master_numParameter, no_dist_numTracesGlobal);

    /* Allocate coordinates on master */
    coordinates.allocate(dist_master_numTracesGlobal);

    /* Local operations on master: 1. Transpose acquisition, 2. calculate 1-D coordinates  */
    if (dist_wavefield->getCommunicator().getRank() == 0) {

        /* Get WriteAccess to local data of acquisition */
        auto write_acquisition_HA = hostWriteAccess(acquisition.getLocalStorage().getData());
        /* Get Readaccess to local data of acquisition_temp */
        auto read_acquisition_temp_HA = hostReadAccess(acquisition_temp.getLocalStorage().getData());

        /* Transpose local data */
        for (IndexType row = 0; row < nrow_temp; row++) {
            for (IndexType column = 0; column < ncolumn_temp; column++) {
                write_acquisition_HA[row + nrow_temp * column] = read_acquisition_temp_HA[column + ncolumn_temp * row];
            }
        }

        /* Release write and read access to local data */
        write_acquisition_HA.release();
        read_acquisition_temp_HA.release();

        /* Get readAccess to acquisition matrix (local) */
        auto read_acquisition_HA = hostReadAccess(acquisition.getLocalStorage().getData());
        /* Get writeAccess to coordinates vector (local) */
        auto write_coordinates_LA = hostWriteAccess(coordinates.getLocalValues());

        Coordinates coord(NX, NY, NZ);

        /* 2. Calculate 1-D coordinates from 3-D coordinates */
        IndexType X, Y, Z;
        for (IndexType i = 0; i < numTracesGlobal; i++) {

            X = read_acquisition_HA[i + numTracesGlobal * 0];
            Y = read_acquisition_HA[i + numTracesGlobal * 1];
            Z = read_acquisition_HA[i + numTracesGlobal * 2];

            write_coordinates_LA[i] = coord.coordinate2index(X, Y, Z);
        }
    }

    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numTracesGlobal);

    /* Get local traces from global traces */
    dmemo::DistributionPtr dist_wavefield_traces = calcDistribution(coordinates, dist_wavefield);

    numTracesLocal = dist_wavefield_traces->getLocalSize();
    numTracesGlobal = dist_wavefield_traces->getGlobalSize();

    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_numParameter, no_dist_numTracesGlobal);

    lama::DenseVector<ValueType> seismogramTypesTmp; // temp vector needed due to type conversion
    acquisition.getRow(seismogramTypesTmp, 3);
    seismogramTypes = lama::cast<IndexType>(seismogramTypesTmp);

    coordinates.redistribute(dist_wavefield_traces);
    seismogramTypes.redistribute(dist_wavefield_traces);

    coordinates.setContextPtr(ctx);
    seismogramTypes.setContextPtr(ctx);

    initOptionalAcquisitionParameter(numParameter, numTracesGlobal, acquisition, dist_wavefield_traces, ctx);
}

/*! \brief Getter method for distribution of local traces 
 *
 \param coordinates coordiantes
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::AcquisitionGeometry<ValueType>::calcDistribution(scai::lama::DenseVector<IndexType> const &coordinates, scai::dmemo::DistributionPtr const dist_wavefield) const
{
    SCAI_ASSERT_DEBUG(coordinates.size() > 0, " The vector coordinates does not contain any elements ! ");

    hmemo::HArray<IndexType> localIndices;

    Global2Local(coordinates, localIndices, dist_wavefield);

    dmemo::DistributionPtr dist_temp(new dmemo::GeneralDistribution(coordinates.size(), localIndices, dist_wavefield->getCommunicatorPtr()));

    return (dist_temp);
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
lama::DenseVector<IndexType> const &KITGPI::Acquisition::AcquisitionGeometry<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_DEBUG(numTracesGlobal == coordinates.size(), "Size mismatch");
    return (coordinates);
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
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::writeAcquisitionToFile(std::string const &filename) const
{
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.assignTranspose(acquisition);
    acquisition_temp.writeToFile(filename);
}

template class KITGPI::Acquisition::AcquisitionGeometry<double>;
template class KITGPI::Acquisition::AcquisitionGeometry<float>;
