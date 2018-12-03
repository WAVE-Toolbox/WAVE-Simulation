
#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "../Configuration/Configuration.hpp"
#include "../Filter/Filter.hpp"
#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include "math.h"
#include "segy.hpp"
#include "suHandler.hpp"

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
            Seismogram() : numSamples(0), numTracesGlobal(0), numTracesLocal(0), normalizeTraces(0), DT(0.0), type(KITGPI::Acquisition::SeismogramType::P){};

            //! Default destructor
            ~Seismogram(){};

            //! Copy Constructor.
            Seismogram(const Seismogram &rhs);

            void swap(KITGPI::Acquisition::Seismogram<ValueType> &rhs);
            void write(Configuration::Configuration const &config, std::string const &filename) const;
            void read(Configuration::Configuration const &config, std::string const &filename, bool copyDist = 0);   
            void read(Configuration::Configuration const &config, std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);
            
            void allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distSeismogram, scai::IndexType NT);

            void redistribute(scai::dmemo::DistributionPtr distRow, scai::dmemo::DistributionPtr distColumn = scai::dmemo::DistributionPtr());
            void replicate();

            void resetData();
            void resetSeismogram();

            void normalizeTrace();
            void integrateTraces();
            void filterTraces(Filter::Filter<ValueType> const &freqFilter);

            /* Getter functions */
            scai::IndexType getNumTracesGlobal() const;
            scai::IndexType getNumTracesLocal() const;
            scai::IndexType getNumSamples() const;
            scai::IndexType getNormalizeTraces() const;
            ValueType getDT() const;
            scai::lama::DenseMatrix<ValueType> &getData();
            scai::lama::DenseMatrix<ValueType> const &getData() const;
            scai::lama::DenseVector<scai::IndexType> const &getCoordinates() const;
            SeismogramType getTraceType() const;
            scai::IndexType getSourceCoordinate() const;
            

            /* Setter functions */
            void setDT(ValueType newDT);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setSourceCoordinate(scai::IndexType sourceCoord);
            void setTraceType(SeismogramType trace);
            void setCoordinates(scai::lama::DenseVector<scai::IndexType> const &coord);
            void setNormalizeTraces(scai::IndexType normalizeTrace);
            void setResampleCoeff(ValueType resampleCoeff = 1.0);

            /* Overloading Operators */
            KITGPI::Acquisition::Seismogram<ValueType> operator*=(ValueType const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator+(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const;
            KITGPI::Acquisition::Seismogram<ValueType> operator+=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> operator-(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const;
            KITGPI::Acquisition::Seismogram<ValueType> operator-=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);
            KITGPI::Acquisition::Seismogram<ValueType> &operator=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs);

          private:
            std::string addSeismogramTypeToName(std::string const &filename) const;
            
            void writeToFileRaw(std::string const &filename) const;
            void writeToFileSU(std::string const &filename, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH) const;
            
            void readFromFileRaw(std::string const &filename, bool copyDist = 0);
            void readFromFileSU(std::string const &filename, bool copyDist = 0);
            
            void readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);
            void readFromFileSU(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);

            scai::IndexType numSamples;      //!< Number of samples of one trace
            scai::IndexType numTracesGlobal; //!< Number of global traces
            scai::IndexType numTracesLocal;  //!< Number of local traces
            scai::IndexType normalizeTraces; //!< L2 Norm of seismogram is calculated

            /* header information */
            ValueType DT;                                         //!< Temporal sampling interval in seconds
            SeismogramType type;                                  //!< Type of trace as #SeismogramType
            scai::lama::DenseVector<scai::IndexType> coordinates; //!< Coordinates of the traces
            scai::IndexType sourceCoordinate;                     //!< Coordinate of source point (in case a single source is used)

            /* raw data */
            scai::lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
            
            /* resampling */
            scai::lama::CSRSparseMatrix<ValueType> resampleMatLeft;
            scai::lama::CSRSparseMatrix<ValueType> resampleMatRight;
            scai::lama::DenseVector<ValueType> resampleVec;
            
        };
    }
}
