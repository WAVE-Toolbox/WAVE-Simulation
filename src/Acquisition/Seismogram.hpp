
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
            Seismogram() : numSamples(0), numTracesGlobal(0), numTracesLocal(0), DT(0.0), type(KITGPI::Acquisition::SeismogramType::P){};

            //! Default destructor
            ~Seismogram(){};

            void write(Configuration::Configuration const &config) const;
            void writeToFileRaw(std::string const &filename) const;
            void writeToFileSU(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, ValueType DH) const;

            void readFromFileRaw(std::string const &filename, dmemo::DistributionPtr distTraces, dmemo::DistributionPtr distSamples);

            void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT);

            void redistribute(dmemo::DistributionPtr distRow, dmemo::DistributionPtr distColumn = NULL);
            void replicate();

            void resetData();

            /* Getter functions */
            IndexType getNumTracesGlobal() const;
            IndexType getNumTracesLocal() const;
            IndexType getNumSamples() const;
            ValueType getDT() const;
            lama::DenseMatrix<ValueType> &getData();
            lama::DenseMatrix<ValueType> const &getData() const;
            lama::DenseVector<IndexType> const &getCoordinates() const;
            SeismogramType getTraceType() const;

            /* Setter functions */
            void setDT(ValueType newDT);
            void setContextPtr(hmemo::ContextPtr ctx);
            void setSourceCoordinate(IndexType sourceCoord);
            void setTraceType(SeismogramType trace);
            void setCoordinates(lama::DenseVector<IndexType> const &coord);

          private:
            std::string addSeismogramTypeToName(std::string const &filename) const;

            IndexType numSamples;      //!< Number of samples of one trace
            IndexType numTracesGlobal; //!< Number of global traces
            IndexType numTracesLocal;  //!< Number of local traces

            /* header information */
            ValueType DT;                             //!< Temporal sampling interval in seconds
            SeismogramType type;                      //!< Type of trace as #SeismogramType
            lama::DenseVector<IndexType> coordinates; //!< Coordinates of the traces
            IndexType sourceCoordinate;               //!< Coordinate of source point (in case a single source is used)

            /* raw data */
            lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
        };
    }
}
