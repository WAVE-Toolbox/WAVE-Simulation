
#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "../Configuration/Configuration.hpp"
#include "../Filter/Filter.hpp"
#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include "math.h"

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
            Seismogram() : DT(0.0), seismoType(KITGPI::Acquisition::SeismogramType::P){};

            //! Default destructor
            ~Seismogram(){};

            //! Copy Constructor.
            Seismogram(const Seismogram &rhs);

            void swap(KITGPI::Acquisition::Seismogram<ValueType> &rhs);
            void write(scai::IndexType const seismogramFormat, std::string const &filename, Coordinates<ValueType> const &modelCoordinates) const;
            void read(scai::IndexType const seismogramFormat, std::string const &filename, bool copyDist = 0);
            //void read(scai::IndexType const SeismogramFormat, std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);

            void allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distSeismogram, scai::IndexType NT);

            void redistribute(scai::dmemo::DistributionPtr distRow, scai::dmemo::DistributionPtr distColumn = scai::dmemo::DistributionPtr());
            void replicate();

            void resetData();
            void resetSeismogram();

            void normalizeTrace();
            void normalizeTraceL2();
            scai::lama::DenseVector<ValueType> getTraceL2norm();
            scai::lama::DenseVector<ValueType> getTraceSum();
            void integrateTraces();
            void filterTraces(Filter::Filter<ValueType> const &freqFilter);

            bool isFinite();
            
            /* Getter functions */
            scai::IndexType getNumTracesGlobal() const;
            scai::IndexType getNumTracesLocal() const;
            scai::IndexType getNumSamples() const;
            ValueType getDT() const;
            scai::lama::DenseMatrix<ValueType> &getData();
            scai::lama::DenseMatrix<ValueType> const &getData() const;
            scai::lama::DenseVector<scai::IndexType> const &get1DCoordinates() const;
            SeismogramType getTraceType() const;
            scai::IndexType getSourceCoordinate() const;

            /* Setter functions */
            void setDT(ValueType newDT);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setSourceCoordinate(scai::IndexType sourceIdx);
            void setTraceType(SeismogramType trace);
            void setCoordinates(scai::lama::DenseVector<scai::IndexType> const &indeces);

            void setSeismoDT(ValueType seismoDT);
            void setEnvelopTrace(scai::IndexType envelopTraces);

            /* Overloading Operators */
            KITGPI::Acquisition::Seismogram<ValueType> operator*=(ValueType const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator*=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator+(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const;
            KITGPI::Acquisition::Seismogram<ValueType> operator+=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator-(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const;
            KITGPI::Acquisition::Seismogram<ValueType> operator-=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> &operator=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);

          private:

            /* header information */
            ValueType DT;                                           //!< Temporal sampling interval in seconds
            ValueType outputDT;                                     //!< Temporal sampling interval for the reampled output in seconds
            SeismogramType seismoType;                              //!< Type of trace as #SeismogramType
            scai::lama::DenseVector<scai::IndexType> coordinates1D; //!< model indeces of the Coordinates of the traces
            scai::IndexType sourceIndex;                            //!< model Index of source point (in case a single source is used)

            /* raw data */
            scai::lama::DenseMatrix<ValueType> data; //!< Raw seismogram data

            /* resampling */
            scai::lama::CSRSparseMatrix<ValueType> resampleMat;
            scai::IndexType outputEnvelope; // output envelope
        };
    }
}
