#pragma once

#include <scai/lama.hpp>

#include "SourceReceiverImplSeismic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace SourceReceiverImpl
        {

            //! \brief FDTDacoustic class
            template <typename ValueType>
            class FDTDacoustic : public SourceReceiverImplSeismic<ValueType>
            {
              public:
                //! Default constructor
                FDTDacoustic(){};
                //! Default destructor
                ~FDTDacoustic(){};

                void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
