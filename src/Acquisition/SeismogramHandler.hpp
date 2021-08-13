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
         * This class stores a Seismogram for each #SeismogramType. This allows easy Seismogram handling, as only one SeismogramHandler can handle a whole multi-component seismogram.
         */
        template <typename ValueType>
        class SeismogramHandler
        {

          public:
            explicit SeismogramHandler();
            //! \brief Default destructor
            ~SeismogramHandler(){};

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
            void setInverseAGC(SeismogramHandler<ValueType> seismograms);

            Seismogram<ValueType> const &getSeismogram(SeismogramType type) const;
            Seismogram<ValueType> &getSeismogram(SeismogramType type);
            scai::IndexType getNumTracesGlobal(SeismogramType type) const;
            Seismogram<ValueType> const &getSeismogram(SeismogramTypeEM type) const;
            Seismogram<ValueType> &getSeismogram(SeismogramTypeEM type);
            scai::IndexType getNumTracesGlobal(SeismogramTypeEM type) const;
            scai::IndexType getNumTracesTotal() const;
            scai::IndexType getNumSamples(SeismogramType type) const;
            scai::IndexType getNumSamples(SeismogramTypeEM type) const;
            KITGPI::Acquisition::SeismogramHandler<ValueType> operator*=(scai::lama::DenseVector<ValueType> const &rhs);
            bool isFinite();

            bool isSeismic = false;
            
          private:
            void setTraceType();
            void setTraceTypeEM();

            std::vector<Seismogram<ValueType>> seismo; //!< Vector which holds the individual Seismogram classes
        };
    }
}
