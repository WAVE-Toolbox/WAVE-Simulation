#pragma once

#include <scai/lama.hpp>

#include "SourceReceiverImpl.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief FDTD2Delastic class
            template <typename ValueType>
            class FDTD2Dsh : public SourceReceiverImpl<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Dsh() = delete;
                //! Default destructor
                ~FDTD2Dsh(){};

                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;

                void applySourcePressure(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("Pressure sources can not be implemented in SH modeling")};
                void applySourceVX(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VX sources can not be implemented in SH modeling")};
                void applySourceVY(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY sources can not be implemented in SH modeling")};
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("Pressure receivers can not be implementedin SH modeling")};
                void gatherSeismogramVX(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VX receivers can not be implementedin SH modeling")};
                void gatherSeismogramVY(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY receivers can not be implementedin SH modeling")};

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
