#pragma once

#include <scai/lama.hpp>

#include "SourceReceiverImplSeismic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief FDTD2Delastic class
            template <typename ValueType>
            class FDTD3Delastic : public SourceReceiverImplSeismic<ValueType>
            {
              public:
                //! Default constructor
                FDTD3Delastic() = delete;
                //! Default destructor
                ~FDTD3Delastic(){};

                using SourceReceiverImplSeismic<ValueType>::SourceReceiverImplSeismic;
                
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
