#pragma once

#include "../Configuration/Configuration.hpp"
#include "Acquisition.hpp"
#include "Seismogram.hpp"
#include <scai/hmemo.hpp>

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

            void readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples);
            void writeToFileRaw(std::string const &filename) const;
            void write(Configuration::Configuration const &config, std::string const &filename) const;
            void normalize();
            void integrate();
            void resetData();

            void setSourceCoordinate(IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setNormalizeTraces(IndexType normalize);
            void setContextPtr(scai::hmemo::ContextPtr ctx);

            Seismogram<ValueType> const &getSeismogram(SeismogramType type) const;
            Seismogram<ValueType> &getSeismogram(SeismogramType type);
            IndexType getNumTracesGlobal(SeismogramType type) const;
            IndexType getNumTracesTotal() const;
            IndexType getNumSamples(SeismogramType type) const;

          private:
            void setTraceType();

            std::vector<Seismogram<ValueType>> seismo; //!< Vector which holds the individual Seismogram classes
        };
    }
}
