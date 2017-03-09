
#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "../Configuration/Configuration.hpp"
#include "Acquisition.hpp"
#include "Coordinates.hpp"

#include "math.h"
#include "segy.h"

namespace KITGPI
{

    namespace Acquisition
    {

        //! Handling of a seismic seismogram
        /*!
         * This class handels a single seismogram which consists of several traces.
         */
        template <typename ValueType>
        class Seismogram
        {

          public:
            //! Default constructor
            Seismogram() : numSamples(0), numTracesGlobal(0), numTracesLocal(0), DT(0.0), normalizeTraces(0), type(KITGPI::Acquisition::SeismogramType::P){};

            //! Default destructor
            ~Seismogram(){};

            //! Copy Constructor.
            Seismogram(const Seismogram &rhs);

            void write(Configuration::Configuration const &config) const;
            void writeToFileRaw(std::string const &filename) const;
            void writeToFileSU(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, ValueType DH) const;

            void readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);

            void allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distSeismogram, IndexType NT);

            void redistribute(scai::dmemo::DistributionPtr distRow, scai::dmemo::DistributionPtr distColumn = NULL);
            void replicate();

            void resetData();

            void normalizeTrace();

            /* Getter functions */
            IndexType getNumTracesGlobal() const;
            IndexType getNumTracesLocal() const;
            IndexType getNumSamples() const;
            IndexType getNormalizeTraces() const;
            ValueType getDT() const;
            scai::lama::DenseMatrix<ValueType> &getData();
            scai::lama::DenseMatrix<ValueType> const &getData() const;
            scai::lama::DenseVector<IndexType> const &getCoordinates() const;
            SeismogramType getTraceType() const;

            /* Setter functions */
            void setDT(ValueType newDT);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setSourceCoordinate(IndexType sourceCoord);
            void setTraceType(SeismogramType trace);
            void setCoordinates(scai::lama::DenseVector<IndexType> const &coord);
            void setNormalizeTraces(IndexType normalizeTrace);

            /* Overloading Operators */
            KITGPI::Acquisition::Seismogram<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator*=(scai::lama::Scalar rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator+(KITGPI::Acquisition::Seismogram<ValueType> rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator+=(KITGPI::Acquisition::Seismogram<ValueType> rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator-(KITGPI::Acquisition::Seismogram<ValueType> rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator-=(KITGPI::Acquisition::Seismogram<ValueType> rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator=(const KITGPI::Acquisition::Seismogram<ValueType> rhs);

          private:
            std::string addSeismogramTypeToName(std::string const &filename) const;

            IndexType numSamples;      //!< Number of samples of one trace
            IndexType numTracesGlobal; //!< Number of global traces
            IndexType numTracesLocal;  //!< Number of local traces

            /* header information */
            ValueType DT;                                   //!< Temporal sampling interval in seconds
            IndexType normalizeTraces;                      //!< L2 Norm of seismogram is calculated
            SeismogramType type;                            //!< Type of trace as #SeismogramType
            scai::lama::DenseVector<IndexType> coordinates; //!< Coordinates of the traces
            IndexType sourceCoordinate;                     //!< Coordinate of source point (in case a single source is used)

            /* raw data */
            scai::lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
        };
    }
}
