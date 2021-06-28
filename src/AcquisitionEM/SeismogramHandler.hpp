#pragma once

#include "../Configuration/Configuration.hpp"
#include "../Filter/Filter.hpp"
#include "Acquisition.hpp"
#include "Seismogram.hpp"
#include <scai/hmemo.hpp>

namespace KITGPI
{

    namespace Acquisition
    {

        /*! \brief Handling of a set of seismograms e.g. different components.
         *
         * This class stores a Seismogram for each #SeismogramTypeEM. This allows easy Seismogram handling, as only one SeismogramHandler can handle a whole multi-component seismogram.
         */
        template <typename ValueType>
        class SeismogramHandlerEM
        {

          public:
            explicit SeismogramHandlerEM();
            //! \brief Default destructor
            ~SeismogramHandlerEM(){};

            void read(scai::IndexType const seismogramFormat, std::string const &filename, bool copyDist = 0);
            void write(scai::IndexType const seismogramFormat, std::string const &filename, Acquisition::Coordinates<ValueType> const &modelCoordinates) const;
            void normalize(scai::IndexType normalizeTraces);
            void integrate();
            void differentiate();
            void resetData();
            void resetSeismograms();
            void filter(Filter::Filter<ValueType> const &freqFilter);

            void setSourceCoordinate(scai::IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setContextPtr(scai::hmemo::ContextPtr ctx);
            void setSeismoDT(ValueType seismoDT);
            void setEnvelopTrace(scai::IndexType envelopTraces);

            SeismogramEM<ValueType> const &getSeismogram(SeismogramTypeEM type) const;
            SeismogramEM<ValueType> &getSeismogram(SeismogramTypeEM type);
            scai::IndexType getNumTracesGlobal(SeismogramTypeEM type) const;
            scai::IndexType getNumTracesTotal() const;
            scai::IndexType getNumSamples(SeismogramTypeEM type) const;
            bool isFinite();

          private:
            void setTraceType();

            std::vector<SeismogramEM<ValueType>> seismo; //!< Vector which holds the individual SeismogramEM classes
        };
    }
}
