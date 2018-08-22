#pragma once

#include "../Configuration/Configuration.hpp"
#include "Acquisition.hpp"
#include "Seismogram.hpp"
#include <scai/hmemo.hpp>
#include "../Filter/Filter.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        /*! \brief Handling of a set of seismograms e.g. different components.
         *
         * This class stores a Seismogram for each #SeismogramType. This allows easy Seismogram handling, as only one SeismogramHandler can handle a whole multi-component seismogram.
         */
        template <typename ValueType>
        class SeismogramHandler
        {

          public:
            explicit SeismogramHandler();
            //! \brief Default destructor
            ~SeismogramHandler(){};

            void readFromFileRaw(std::string const &filename, bool copyDist = 0);
            void readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);
            void writeToFileRaw(std::string const &filename) const;
            void write(Configuration::Configuration const &config, std::string const &filename) const;
            void normalize();
            void integrate();
            void resetData();
            void filter(Filter::Filter<ValueType> const &freqFilter);

            void setSourceCoordinate(scai::IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setNormalizeTraces(scai::IndexType normalize);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setResampleCoeff(scai::IndexType resampleCoeff = 1);

            Seismogram<ValueType> const &getSeismogram(SeismogramType type) const;
            Seismogram<ValueType> &getSeismogram(SeismogramType type);
            scai::IndexType getNumTracesGlobal(SeismogramType type) const;
            scai::IndexType getNumTracesTotal() const;
            scai::IndexType getNumSamples(SeismogramType type) const;

          private:
            void setTraceType();

            std::vector<Seismogram<ValueType>> seismo; //!< Vector which holds the individual Seismogram classes
        };
    }
}
