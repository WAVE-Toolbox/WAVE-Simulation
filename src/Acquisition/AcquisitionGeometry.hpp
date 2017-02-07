#pragma once

#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Acquisition.hpp"

#include "Coordinates.hpp"

#include "SeismogramHandler.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        template <typename ValueType>
        class AcquisitionGeometry
        {

          public:
            AcquisitionGeometry() : numTracesGlobal(0), numTracesLocal(0), numParameter(0){};
            ~AcquisitionGeometry(){};

            /* I/O for acquisition */
            void readAcquisitionFromFile(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx);
            void writeAcquisitionToFile(std::string const &filename) const;

            /* Getter methods */
            lama::DenseVector<IndexType> const &getCoordinates() const;
            lama::DenseVector<IndexType> const &getSeismogramTypes() const;
            IndexType getNumTracesGlobal() const;
            IndexType getNumTracesLocal() const;
            SeismogramHandler<ValueType> &getSeismogramHandler();
            SeismogramHandler<ValueType> const &getSeismogramHandler() const;

          protected:
            IndexType getNumParameter() const;
            void initSeismogramHandler(IndexType const NT, hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield);

          private:
            IndexType numTracesGlobal; //!< Number of global traces
            IndexType numTracesLocal;  //!< Number of local traces

            /* SeismogramHandler to store seismic data */
            SeismogramHandler<ValueType> seismograms; //!< SeismogramHandler to handle the #Seismogram's

            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition;     //!< Matrix that stores the acquisition
            IndexType numParameter;                       //!< Number of parameters given in acquisition matrix
            lama::DenseVector<IndexType> coordinates;     //!< Global coordinates of the traces (1-D coordinates)
            lama::DenseVector<IndexType> seismogramTypes; //!< #SeismogramType of the traces: 1==Pressure, 2==vX, 3==vY, 4==vZ

            /* Methods in derived classes for the readAcquisitionFromFile method */
            virtual void checkRequiredNumParameter(IndexType numParameterCheck) = 0;
            virtual void initOptionalAcquisitionParameter(IndexType numParameter, IndexType numTracesGlobal, lama::DenseMatrix<ValueType> acquisition, dmemo::DistributionPtr dist_wavefield_traces, hmemo::ContextPtr ctx);

            /* Calculation of distribution for local  traces */
            dmemo::DistributionPtr calcDistribution(lama::DenseVector<IndexType> const &coordinates, dmemo::DistributionPtr const dist_wavefield) const;
        };
    }
}

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initOptionalAcquisitionParameter(IndexType /*numParameter*/, IndexType /*numTracesGlobal*/, lama::DenseMatrix<ValueType> /*acquisition*/, dmemo::DistributionPtr /*dist_wavefield_traces*/, hmemo::ContextPtr /*ctx*/)
{
}

template <typename ValueType>
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::initSeismogramHandler(IndexType const NT, hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield)
{

    SCAI_ASSERT_DEBUG(seismogramTypes.size() == coordinates.size(), "Size mismatch");
    SCAI_ASSERT_DEBUG(numTracesGlobal == seismogramTypes.size(), "Size mismatch");

    IndexType count[NUM_ELEMENTS_SEISMOGRAMTYPE] = {0, 0, 0, 0};
    lama::DenseVector<IndexType> coord[NUM_ELEMENTS_SEISMOGRAMTYPE];
    dmemo::DistributionPtr dist[NUM_ELEMENTS_SEISMOGRAMTYPE];

    /* Count elements for each source type */
    lama::Scalar tempScalar;
    IndexType tempIndexType;
    for (IndexType i = 0; i < numTracesGlobal; ++i) {
        tempScalar = seismogramTypes.getValue(i);
        tempIndexType = tempScalar.getValue<IndexType>() - 1;
        SCAI_ASSERT_DEBUG(tempIndexType >= 0 && tempIndexType <= 3, "Unkown Source Type");

        ++count[tempIndexType];
    }

    /* Allocate lama vectors */
    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; ++i) {
        coord[i].allocate(count[i]);
        count[i] = 0;
    }

    /* Sort coordinates */
    for (IndexType i = 0; i < numTracesGlobal; ++i) {

        tempScalar = seismogramTypes.getValue(i);
        tempIndexType = tempScalar.getValue<IndexType>() - 1;

        coord[tempIndexType].setValue(count[tempIndexType], this->getCoordinates().getValue(i));
        ++count[tempIndexType];
    }

    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(0) == SeismogramType::P, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(1) == SeismogramType::VX, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(2) == SeismogramType::VY, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(3) == SeismogramType::VZ, "Cast went wrong");

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
void KITGPI::Acquisition::AcquisitionGeometry<ValueType>::readAcquisitionFromFile(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx)
{

    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NX<=0");

    /* Read acquisition matrix */
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(filename);

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
        lama::DenseStorage<ValueType> *acquisition_DS = &acquisition.getLocalStorage();
        hmemo::HArray<ValueType> *acquisition_HA = &acquisition_DS->getData();
        hmemo::WriteAccess<ValueType> write_acquisition_HA(*acquisition_HA);

        /* Get Readaccess to local data of acquisition_temp */
        lama::DenseStorage<ValueType> *acquisition_temp_DS = &acquisition_temp.getLocalStorage();
        hmemo::HArray<ValueType> *acquisition_temp_HA = &acquisition_temp_DS->getData();
        hmemo::ReadAccess<ValueType> read_acquisition_temp_HA(*acquisition_temp_HA);

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
        acquisition_DS = &acquisition.getLocalStorage();
        acquisition_HA = &acquisition_DS->getData();
        hmemo::ReadAccess<ValueType> read_acquisition_HA(*acquisition_HA);

        /* Get writeAccess to coordinates vector (local) */
        utilskernel::LArray<IndexType> *coordinates_LA = &coordinates.getLocalValues();
        hmemo::WriteAccess<IndexType> write_coordinates_LA(*coordinates_LA);

        Coordinates<ValueType> coord;

        /* 2. Calculate 1-D coordinates form 3-D coordinates */
        IndexType X, Y, Z;
        for (IndexType i = 0; i < numTracesGlobal; i++) {

            X = read_acquisition_HA[i + numTracesGlobal * 0];
            Y = read_acquisition_HA[i + numTracesGlobal * 1];
            Z = read_acquisition_HA[i + numTracesGlobal * 2];

            write_coordinates_LA[i] = coord.coordinate2index(X, Y, Z, NX, NY, NZ);
        }

        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
    }

    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numTracesGlobal);

    /* Get local traces from global traces */
    dmemo::DistributionPtr dist_wavefield_traces = calcDistribution(coordinates, dist_wavefield);

    numTracesLocal = dist_wavefield_traces->getLocalSize();
    numTracesGlobal = dist_wavefield_traces->getGlobalSize();

    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_numParameter, no_dist_numTracesGlobal);

    seismogramTypes.allocate(numTracesGlobal);
    acquisition.getRow(seismogramTypes, 3);

    coordinates.redistribute(dist_wavefield_traces);
    seismogramTypes.redistribute(dist_wavefield_traces);

    coordinates.setContextPtr(ctx);
    seismogramTypes.setContextPtr(ctx);

    initOptionalAcquisitionParameter(numParameter, numTracesGlobal, acquisition, dist_wavefield_traces, ctx);
}

/*! \brief Getter methode for Distribution.
 *
 \param coordinates coordiantes
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::AcquisitionGeometry<ValueType>::calcDistribution(lama::DenseVector<IndexType> const &coordinates, dmemo::DistributionPtr const dist_wavefield) const
{
    SCAI_ASSERT_DEBUG(coordinates.size() > 0, " The vector coordinates does not contain any elements ! ");

    hmemo::HArray<IndexType> localIndices;

    Coordinates<ValueType> coord;
    coord.Global2Local(coordinates, localIndices, dist_wavefield);

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
