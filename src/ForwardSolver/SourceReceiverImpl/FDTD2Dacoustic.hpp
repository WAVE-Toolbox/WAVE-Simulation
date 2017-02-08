
#pragma once

#include <scai/lama.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

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

                using FDTDacoustic<ValueType>::FDTDacoustic;
            };
        }
    }
}
