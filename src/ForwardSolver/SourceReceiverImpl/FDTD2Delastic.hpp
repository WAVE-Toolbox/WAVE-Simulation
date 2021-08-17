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
            class FDTD2Delastic : public SourceReceiverImplSeismic<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Delastic() = delete;
                //! Default destructor
                ~FDTD2Delastic(){};

                using SourceReceiverImplSeismic<ValueType>::SourceReceiverImplSeismic;
                
                void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefieldIN, scai::IndexType t) override;

                void applySourceVZ(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY sources can not be implemented in 2D elastic modeling")};
                void gatherSeismogramVZ(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY receivers can not be implemented in 2D elastic modeling")};

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
