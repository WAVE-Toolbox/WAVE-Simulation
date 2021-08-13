#pragma once

#include "../Configuration/Configuration.hpp"
#include "../Filter/Filter.hpp"
#include "../Acquisition/Acquisition.hpp"
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
            void write(scai::IndexType const seismogramFormat, std::string const &filename, Coordinates<ValueType> const &modelCoordinates) const;
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
            void setInstantaneousTrace(scai::IndexType instantaneousTraces);
            void setFrequencyAGC(ValueType setFrequencyAGC);
            void calcInverseAGC();
            void setInverseAGC(SeismogramHandlerEM<ValueType> seismograms);

            SeismogramEM<ValueType> const &getSeismogram(SeismogramTypeEM type) const;
            SeismogramEM<ValueType> &getSeismogram(SeismogramTypeEM type);
            scai::IndexType getNumTracesGlobal(SeismogramTypeEM type) const;
            scai::IndexType getNumTracesTotal() const;
            scai::IndexType getNumSamples(SeismogramTypeEM type) const;
            KITGPI::Acquisition::SeismogramHandlerEM<ValueType> operator*=(scai::lama::DenseVector<ValueType> const &rhs);
            bool isFinite();

          private:
            void setTraceTypeEM();

            std::vector<SeismogramEM<ValueType>> seismo; //!< Vector which holds the individual SeismogramEM classes
        };
    }
}
