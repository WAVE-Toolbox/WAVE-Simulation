#pragma once

#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Acquisition.hpp"
#include "Coordinates.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        template <typename ValueType>
        class SeismogramHandler;

        /*! \brief Handling of receivers
         *
         * This class accounts for the handling of seismic receivers.
         * It provides the reading of the receivers acquisition from file, the distribution of the receivers and the collection of the seismograms.
         *
         *
         * This class will first read-in the global receiver configuration from a receiver configuration file and
         * afterwards it will determine the individual #SeismogramType of each receiver, and based on the #SeismogramType it will
         * initialize the #Seismogram's of the #SeismogramHandler #receiver.\n\n
         * The #coordinates and #receiver_type vectors contain the coordinates and #SeismogramType of all receivers, whereas the
         * individual Seismogram of the SeismogramHandler #receiver will hold the receivers seperated based on the #SeismogramType.
         *
         */
        template <typename ValueType>
        class Receivers
        {
          public:
            //! \brief Default constructor
            Receivers() : numReceiversGlobal(0), numReceiversLocal(0), numParameter(0){};
            explicit Receivers(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);

            //! \brief Default destructor
            ~Receivers(){};

            void init(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);

            void readReceiverAcquisition(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx);
            void writeReceiverAcquisition(std::string const &filename) const;

            /* Getter functions */
            IndexType getNumReceiversGlobal() const;
            IndexType getNumReceiversLocal() const;
            dmemo::DistributionPtr getReceiversDistribution() const;
            lama::DenseVector<IndexType> const &getCoordinates() const;
            lama::DenseVector<IndexType> const &getReceiversType() const;
            SeismogramHandler<ValueType> &getSeismogramHandler();

          private:
            void initSeismogramHandler(IndexType const NT, hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield);

            dmemo::DistributionPtr getReceiverDistribution(lama::DenseVector<IndexType> const &coordinates, dmemo::DistributionPtr const dist_wavefield) const;

            IndexType numReceiversGlobal; //!< Number of receivers global
            IndexType numReceiversLocal;  //!< Number of receivers local

            SeismogramHandler<ValueType> receiver; //!< SeismogramHandler to handle the receiver seismograms

            dmemo::DistributionPtr dist_wavefield_receivers; //!< Calculated Distribution of the receivers based on the distribution of the wavefields
            dmemo::DistributionPtr no_dist_NT;               //!< No distribution of the columns of the seismogram matrix

            /* Acquisition Settings */
            lama::DenseMatrix<ValueType> acquisition;   //!< Matrix that stores the receiver acquisition
            IndexType numParameter;                     //!< Number of receiver parameters given in acquisition matrix
            lama::DenseVector<IndexType> coordinates;   //!< Coordinates of receivers global (1-D coordinates)
            lama::DenseVector<IndexType> receiver_type; //!< #SeismogramType of the receivers: 1==Pressure, 2==vX, 3==vY, 4==vZ
        };
    }
}

/*! \brief Initialization of the receiver SeismogramHandler
 *
 * This method initializes the SeismogramHandler, which is used to save the seismograms of the receivers.
 * This method will determine which receivers belong to which #SeismogramType. 
 * Based on the receiveres for each #SeismogramType the method will calculate the distribution for the receiveres for each #SeismogramType and will use this distribution to allocate the #SeismogramHandler.
 \param NT Numer of timesteps
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::initSeismogramHandler(IndexType const NT, hmemo::ContextPtr const ctx, dmemo::DistributionPtr const dist_wavefield)
{

    SCAI_ASSERT_DEBUG(receiver_type.size() == coordinates.size(), "Size mismatch")

    IndexType count[NUM_ELEMENTS_SEISMOGRAMTYPE] = {0, 0, 0, 0};
    lama::DenseVector<IndexType> coord[NUM_ELEMENTS_SEISMOGRAMTYPE];
    dmemo::DistributionPtr dist[NUM_ELEMENTS_SEISMOGRAMTYPE];

    IndexType numReceiversGlobal = receiver_type.size();

    /* Count elements for each source type */
    lama::Scalar tempScalar;
    IndexType tempIndexType;
    for (IndexType i = 0; i < numReceiversGlobal; ++i) {
        tempScalar = receiver_type.getValue(i);
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
    for (IndexType i = 0; i < numReceiversGlobal; ++i) {

        tempScalar = receiver_type.getValue(i);
        tempIndexType = tempScalar.getValue<IndexType>() - 1;

        coord[tempIndexType].setValue(count[tempIndexType], coordinates.getValue(i));
        ++count[tempIndexType];
    }

    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(0) == SeismogramType::P, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(1) == SeismogramType::VX, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(2) == SeismogramType::VY, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(3) == SeismogramType::VZ, "Cast went wrong");

    /* Calculate distribution, redistribute coordinates and set coordinates to seismogramHandler */
    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; ++i) {

        if (coord[i].size() > 0) {
            dist[i] = getReceiverDistribution(coord[i], dist_wavefield);
            receiver.getSeismogram(static_cast<SeismogramType>(i)).allocate(ctx, dist[i], NT);
            coord[i].redistribute(dist[i]);
            receiver.getSeismogram(static_cast<SeismogramType>(i)).setCoordinates(coord[i]);
        }
        count[i] = 0;
    }

    receiver.setContextPtr(ctx);
}

/*! \brief Getter method for seismogram handler
 * 
 * This method will return the SeismogramHandler of this receiver class. 
 \return The SeismogramHandler of the receivers.
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType> &KITGPI::Acquisition::Receivers<ValueType>::getSeismogramHandler()
{
    return (receiver);
}

/*! \brief Getter method for reference to the vector containing the receiver type
 *
 * The values of the receiver types correspond to the #SeismogramType enum.\n
 * The size of the receiver type vector is equal to the number of global receivers.
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::Receivers<ValueType>::getReceiversType() const
{
    SCAI_ASSERT_ERROR(receiver_type.size() != 0, "No receivers type set ");
    return (receiver_type);
}

/*! \brief Getter method for reference to the coordinates vector
 * 
 * The coordinates vector stores the 1-D coordinates of the receivers.\n
 * The size of the coordinates vector is equal to the number of global receivers.
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::Receivers<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_ERROR(coordinates.size() != 0, "No receivers coordinates set ");
    return (coordinates);
}

/*! \brief Getter method for the global distribution of the receivers based on the wavefield distribution
 *
 * This method returns the distribution of the global receivers.
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::Receivers<ValueType>::getReceiversDistribution() const
{
    SCAI_ASSERT_ERROR(dist_wavefield_receivers != nullptr, "Receivers distribution not set ");
    return (dist_wavefield_receivers);
}

/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
KITGPI::Acquisition::Receivers<ValueType>::Receivers(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
    : numReceiversGlobal(0), numReceiversLocal(0), numParameter(0)
{
    init(config, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
{
    readReceiverAcquisition(config.get<std::string>("ReceiverFilename"), config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), dist_wavefield, ctx);
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
    initSeismogramHandler(getNT, ctx, dist_wavefield);
    receiver.setDT(config.get<ValueType>("DT"));
}

/*! \brief Get number of global receivers
 *
 \return Number of global receivers
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversGlobal() const
{
    return (numReceiversGlobal);
}

/*! \brief Get number of local receivers
 *
 \return Number of local receivers on this process
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Receivers<ValueType>::getNumReceiversLocal() const
{
    return (numReceiversLocal);
}

/*! \brief Write receivers acquisition to file
 *
 \param filename Filename to write receivers acquisition
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::writeReceiverAcquisition(std::string const &filename) const
{
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.assignTranspose(acquisition);
    acquisition_temp.writeToFile(filename);
}

/*! \brief Read receivers acquisition from file
 *
 * This method reads the receivers acquisition, calculates the 1-D coordinates from the 3-D coordinates,
 * splits up the receivers configuration into the corresponding vectors and calculates the receivers distribution.
 * The parameter vectors will be distributed accordingly to the receiver distribution.
 *
 \param filename Filename to read receivers acquisition
 \param NX Number of global grid points in X
 \param NY Number of global grid points in Y
 \param NZ Number of global grid points in Z
 \param dist_wavefield Distribution of the wavefields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::readReceiverAcquisition(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist_wavefield, hmemo::ContextPtr ctx)
{

    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NX<=0");

    /* Read acquisition matrix */
    lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(filename);

    IndexType nrow_temp = acquisition_temp.getNumRows();
    IndexType ncolumn_temp = acquisition_temp.getNumColumns();

    /* Derive number of receivers and number of read-in parameters */
    numReceiversGlobal = nrow_temp;
    numParameter = ncolumn_temp;

    /* Check if number of parameters is supported */
    if (numParameter != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }

    /* Distribution: Master process only (comm->myRank()==0) */
    dmemo::DistributionPtr dist_master_numParameter(new dmemo::CyclicDistribution(numParameter, numParameter, dist_wavefield->getCommunicatorPtr()));
    dmemo::DistributionPtr dist_master_numReceiversGlobal(new dmemo::CyclicDistribution(numReceiversGlobal, numReceiversGlobal, dist_wavefield->getCommunicatorPtr()));

    /* Distribution: Replicated on all processes */
    dmemo::DistributionPtr no_dist_numReceiversGlobal(new scai::dmemo::NoDistribution(numReceiversGlobal));
    dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(numParameter));

    /* Allocate acquisition matrix on master */
    acquisition.allocate(dist_master_numParameter, no_dist_numReceiversGlobal);

    /* Allocate coordinates on master */
    coordinates.allocate(dist_master_numReceiversGlobal);

    /* Allocate receiver parameter vectors on master */
    receiver_type.allocate(dist_master_numReceiversGlobal);

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
        for (IndexType i = 0; i < numReceiversGlobal; i++) {

            X = read_acquisition_HA[i + numReceiversGlobal * 0];
            Y = read_acquisition_HA[i + numReceiversGlobal * 1];
            Z = read_acquisition_HA[i + numReceiversGlobal * 2];

            write_coordinates_LA[i] = coord.coordinate2index(X, Y, Z, NX, NY, NZ);
        }

        /* Release write and read access to local data */
        read_acquisition_HA.release();
        write_coordinates_LA.release();
    }

    /* Replicate coordinates on all processes */
    coordinates.redistribute(no_dist_numReceiversGlobal);

    dist_wavefield_receivers = getReceiverDistribution(coordinates, dist_wavefield);
    numReceiversLocal = dist_wavefield_receivers->getLocalSize();
    numReceiversGlobal = dist_wavefield_receivers->getGlobalSize();

    /* Replicate acquisition on all processes */
    acquisition.redistribute(no_dist_numParameter, no_dist_numReceiversGlobal);

    /* Allocate receiver parameter vectors on all processes */
    receiver_type.allocate(numReceiversGlobal);

    /* Save receiver configurations from acquisition matrix in vectors */
    acquisition.getRow(receiver_type, 3);

    /* Redistribute receiver parameter vectors to corresponding processes */
    coordinates.redistribute(dist_wavefield_receivers);
    receiver_type.redistribute(dist_wavefield_receivers);

    coordinates.setContextPtr(ctx);
    receiver_type.setContextPtr(ctx);
}

/*! \brief Calculation of the distribution of the receivers
 *
 * This method calculates the local receivers based on the vector conatining the global receiver coordinates. 
 \param coordinates Filename to read receivers acquisition
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Acquisition::Receivers<ValueType>::getReceiverDistribution(lama::DenseVector<IndexType> const &coordinates, dmemo::DistributionPtr const dist_wavefield) const
{
    SCAI_ASSERT_DEBUG(coordinates.size() > 0, " The vector coordinates does not contain any elements ! ");

    hmemo::HArray<IndexType> localIndices;

    Coordinates<ValueType> coord;
    coord.Global2Local(coordinates, localIndices, dist_wavefield);

    dmemo::DistributionPtr dist_temp(new dmemo::GeneralDistribution(coordinates.size(), localIndices, dist_wavefield->getCommunicatorPtr()));

    return (dist_temp);
}
