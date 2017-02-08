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

