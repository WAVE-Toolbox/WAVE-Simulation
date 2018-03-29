#pragma once

#include <scai/lama.hpp>

#include "FDTDacoustic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace SourceReceiverImpl
        {

            //! \brief FDTD2Dacoustic class
            template <typename ValueType>
            class FDTD2Dacoustic : public FDTDacoustic<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Dacoustic() = delete;
                //! Default destructor
                ~FDTD2Dacoustic(){};
		
                void applySourceVZ(Acquisition::Seismogram<ValueType> const &/*seismo*/, Wavefields::Wavefields<ValueType> &/*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY sources can not be implemented in 2D acoustic modeling")};
                void gatherSeismogramVZ(Acquisition::Seismogram<ValueType> &/*seismo*/, Wavefields::Wavefields<ValueType> &/*wavefield*/, IndexType /*t*/) override{COMMON_THROWEXCEPTION("VY receivers can not be implementedin 2D acoustic modeling")};

                using FDTDacoustic<ValueType>::FDTDacoustic;
            };
        }
    }
}
