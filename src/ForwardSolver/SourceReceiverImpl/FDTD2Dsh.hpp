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
            class FDTD2Dsh : public SourceReceiverImplSeismic<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Dsh() = delete;
                //! Default destructor
                ~FDTD2Dsh(){};

                using SourceReceiverImplSeismic<ValueType>::SourceReceiverImplSeismic;
                
                void applySourcePressure(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("Pressure sources can not be implemented in SH modeling")};
                void applySourceVX(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VX sources can not be implemented in SH modeling")};
                void applySourceVY(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY sources can not be implemented in SH modeling")};
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("Pressure receivers can not be implemented in SH modeling")};
                void gatherSeismogramVX(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VX receivers can not be implemented in SH modeling")};
                void gatherSeismogramVY(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY receivers can not be implemented in SH modeling")};

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
