#pragma once

#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <vector>

#include "Acquisition.hpp"

#include "Coordinates.hpp"

#include "AcquisitionSettings.hpp"
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
            template <class settingsVec>
            void setAcquisition(settingsVec allSettings, Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx);
            void writeAcquisitionToFile(std::string const &filename) const;

            /* Getter methods */
            scai::lama::DenseVector<scai::IndexType> const &get1DCoordinates() const;
            scai::lama::DenseVector<scai::IndexType> const &getSeismogramTypes() const;
            scai::IndexType getNumTracesGlobal() const;
            scai::IndexType getNumTracesLocal() const;
            SeismogramHandler<ValueType> &getSeismogramHandler();
            SeismogramHandler<ValueType> const &getSeismogramHandler() const;

          protected:
            scai::IndexType getNumParameter() const;
            void initSeismogramHandler(scai::IndexType const NT, scai::hmemo::ContextPtr const ctx, scai::dmemo::DistributionPtr const dist_wavefield);

          private:
            scai::IndexType numTracesGlobal; //!< Number of global traces
            scai::IndexType numTracesLocal;  //!< Number of local traces

            /* SeismogramHandler to store seismic data */
            SeismogramHandler<ValueType> seismograms; //!< SeismogramHandler to handle the #Seismogram's

            /* Acquisition Settings */
            scai::IndexType numParameter;                             //!< Number of parameters given in acquisition matrix
            scai::lama::DenseVector<scai::IndexType> coordinates1D;   //!< Global coordinates of the traces (1-D coordinates)
            scai::lama::DenseVector<scai::IndexType> seismogramTypes; //!< #SeismogramType of the traces: 1==Pressure, 2==vX, 3==vY, 4==vZ

            /* Methods in derived classes for the readAcquisitionFromFile method */
            virtual void checkRequiredNumParameter(scai::IndexType numParameterCheck) = 0;
            virtual void initOptionalAcquisitionParameter(scai::IndexType numParameter, scai::IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx);

        };

        /*! \brief reads parameters from the acquisition matrix and redistributes 
        *
        * seismogram coordinates and Seismogram Types are stored in vectors and redistributet according to the distribution of the wavefield.
        * seismograms are only stored at the spatial domains which include the receiver/source coordinate.
        \param acquisition_matrix Acquisition matrix (contains eg. seismogram coordinates and seismogram types)
        \param modelCoordinates Coordinates object (handles wavefield coordinates)
        \param dist_wavefield Distribution of the wavefield
        \param ctx Context
        */
        template <typename ValueType>
        template <class settingsVec>
        void AcquisitionGeometry<ValueType>::setAcquisition(settingsVec allSettings, Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx)
        {
            scai::IndexType nrow_temp = allSettings.size();

            /* Derive number of traces and number of read-in parameters */
            numTracesGlobal = nrow_temp;

            /* Distribution: Master process only (comm->myRank()==0) */
            scai::dmemo::DistributionPtr dist_master_numTracesGlobal(new scai::dmemo::CyclicDistribution(numTracesGlobal, numTracesGlobal, dist_wavefield->getCommunicatorPtr()));

            /* Distribution: Replicated on all processes */
            scai::dmemo::DistributionPtr no_dist_numTracesGlobal(new scai::dmemo::NoDistribution(numTracesGlobal));

            /* Allocate coordinates on master */
            coordinates1D.allocate(dist_master_numTracesGlobal);
            seismogramTypes.allocate(dist_master_numTracesGlobal);

            /* Local operations on master: 1. Transpose acquisition, 2. calculate 1-D coordinates  */
            if (dist_wavefield->getCommunicator().getRank() == 0) {
                /* Get writeAccess to coordinates vector (local) */
                auto write_coordinates_LA = hostWriteAccess(coordinates1D.getLocalValues());
                auto write_seismogramTypes_LA = hostWriteAccess(seismogramTypes.getLocalValues());

                /* 2. Calculate 1-D coordinates from 3-D coordinates */
                for (scai::IndexType i = 0; i < numTracesGlobal; i++) {
                    write_coordinates_LA[i] = modelCoordinates.coordinate2index(allSettings[i].getCoords());
                    write_seismogramTypes_LA[i] = allSettings[i].getType();
                }
            }

            /* Replicate coordinates on all processes */
            coordinates1D.redistribute(no_dist_numTracesGlobal);

            /* Get local traces from global traces */
            scai::dmemo::DistributionPtr dist_wavefield_traces = calcDistribution(coordinates1D, dist_wavefield);

            numTracesLocal = dist_wavefield_traces->getLocalSize();
            numTracesGlobal = dist_wavefield_traces->getGlobalSize();

            coordinates1D.redistribute(dist_wavefield_traces);
            seismogramTypes.redistribute(dist_wavefield_traces);

            coordinates1D.setContextPtr(ctx);
            seismogramTypes.setContextPtr(ctx);
        }
    }
}
