
#pragma once

#include <scai/lama.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

#include "SourceReceiverImpl.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace SourceReceiverImpl
        {

            //! \brief FDTDacoustic class
            template <typename ValueType>
            class FDTDacoustic : public SourceReceiverImpl<ValueType>
            {
              public:
                //! Default constructor
                FDTDacoustic() = delete;
                //! Default destructor
                ~FDTDacoustic(){};

                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;

                void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t) override;
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t) override;

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
