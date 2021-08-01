
#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "../Configuration/Configuration.hpp"
#include "../Filter/Filter.hpp"
#include "Acquisition.hpp"
#include "../Acquisition/Coordinates.hpp"
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
        class SeismogramEM
        {

          public:
            //! Default constructor
            SeismogramEM() : DT(0.0), seismoType(KITGPI::Acquisition::SeismogramTypeEM::EZ){};

            //! Default destructor
            ~SeismogramEM(){};

            //! Copy Constructor.
            SeismogramEM(const SeismogramEM &rhs);

            void swap(KITGPI::Acquisition::SeismogramEM<ValueType> &rhs);
            void write(scai::IndexType const seismogramFormat, std::string const &filename, Acquisition::Coordinates<ValueType> const &modelCoordinates) const;
            void read(scai::IndexType const seismogramFormat, std::string const &filename, bool copyDist = 0);
            //void read(scai::IndexType const SeismogramFormat, std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);

            void allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distSeismogram, scai::IndexType NT);

            void redistribute(scai::dmemo::DistributionPtr distRow, scai::dmemo::DistributionPtr distColumn = scai::dmemo::DistributionPtr());
            void replicate();

            void resetData();
            void resetSeismogram();

            void normalizeTrace(scai::IndexType normalizeTraces);
            scai::lama::DenseMatrix<ValueType> getAGCSum();
            void calcInverseAGC();
            scai::lama::DenseMatrix<ValueType> const &getInverseAGC() const;
            void setInverseAGC(scai::lama::DenseMatrix<ValueType> setInverseAGC);
            scai::lama::DenseVector<ValueType> getTraceL2norm();
            scai::lama::DenseVector<ValueType> getTraceSum();
            void integrateTraces();
            void differentiateTraces();
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
            SeismogramTypeEM getTraceType() const;
            scai::IndexType getSourceCoordinate() const;

            /* Setter functions */
            void setDT(ValueType newDT);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setSourceCoordinate(scai::IndexType sourceIdx);
            void setTraceType(SeismogramTypeEM trace);
            void setCoordinates(scai::lama::DenseVector<scai::IndexType> const &indeces);

            void setSeismoDT(ValueType seismoDT);
            void setInstantaneousTrace(scai::IndexType instantaneousTraces);
            void setFrequencyAGC(ValueType setFrequencyAGC);

            /* Overloading Operators */
            KITGPI::Acquisition::SeismogramEM<ValueType> operator*=(ValueType const &rhs);
            KITGPI::Acquisition::SeismogramEM<ValueType> operator*=(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs);
            KITGPI::Acquisition::SeismogramEM<ValueType> operator*=(scai::lama::DenseVector<ValueType> const &rhs);
            KITGPI::Acquisition::SeismogramEM<ValueType> operator+(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs) const;
            KITGPI::Acquisition::SeismogramEM<ValueType> operator+=(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs);
            KITGPI::Acquisition::SeismogramEM<ValueType> operator-(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs) const;
            KITGPI::Acquisition::SeismogramEM<ValueType> operator-=(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs);
            KITGPI::Acquisition::SeismogramEM<ValueType> &operator=(KITGPI::Acquisition::SeismogramEM<ValueType> const &rhs);

          private:

            /* header information */
            ValueType DT;                                           //!< Temporal sampling interval in seconds
            ValueType outputDT;                                     //!< Temporal sampling interval for the reampled output in seconds
            SeismogramTypeEM seismoType;                                    //!< Type of trace as #SeismogramTypeEM
            scai::lama::DenseVector<scai::IndexType> coordinates1D; //!< model indeces of the Coordinates of the traces
            scai::IndexType sourceIndex;                            //!< model Index of source point (in case a single source is used)

            /* raw data */
            scai::lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
            scai::lama::DenseMatrix<ValueType> inverseAGC; //!< inverse of AGC

            /* resampling */
            scai::lama::CSRSparseMatrix<ValueType> resampleMat;
            scai::IndexType outputInstantaneous; // output envelope
            ValueType frequencyAGC; // frequency used to calculate AGC window length
            bool useAGC = true; // make sure AGC can be applied only once for each shot
        };
    }
}
