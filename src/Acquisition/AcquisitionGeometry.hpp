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
            void setAcquisition(scai::lama::DenseMatrix<ValueType> acquisition_temp, Coordinates const &modelCoordinates, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx);
            void writeAcquisitionToFile(std::string const &filename) const;

            /* Getter methods */
            scai::lama::DenseVector<scai::IndexType> const &getCoordinates() const;
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

            void Global2Local(scai::lama::Vector<scai::IndexType> const &coordinatesglobal, scai::hmemo::HArray<scai::IndexType> &coordinateslocal, scai::dmemo::DistributionPtr dist) const;

            /* SeismogramHandler to store seismic data */
            SeismogramHandler<ValueType> seismograms; //!< SeismogramHandler to handle the #Seismogram's

            /* Acquisition Settings */
            scai::lama::DenseMatrix<ValueType> acquisition;           //!< Matrix that stores the acquisition
            scai::IndexType numParameter;                             //!< Number of parameters given in acquisition matrix
            scai::lama::DenseVector<scai::IndexType> coordinates;     //!< Global coordinates of the traces (1-D coordinates)
            scai::lama::DenseVector<scai::IndexType> seismogramTypes; //!< #SeismogramType of the traces: 1==Pressure, 2==vX, 3==vY, 4==vZ

            /* Methods in derived classes for the readAcquisitionFromFile method */
            virtual void checkRequiredNumParameter(scai::IndexType numParameterCheck) = 0;
            virtual void initOptionalAcquisitionParameter(scai::IndexType numParameter, scai::IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx);

            /* Calculation of distribution for local  traces */
            scai::dmemo::DistributionPtr calcDistribution(scai::lama::DenseVector<scai::IndexType> const &coordinates, scai::dmemo::DistributionPtr const dist_wavefield) const;
        };
    }
}
