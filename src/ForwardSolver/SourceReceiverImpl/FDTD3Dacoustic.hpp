
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

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

            //! \brief FDTD3Dacoustic class
            template <typename ValueType>
            class FDTD3Dacoustic : public FDTDacoustic<ValueType>
            {
              public:
                //! Default constructor
                FDTD3Dacoustic() = delete;
                //! Default destructor
                ~FDTD3Dacoustic(){};

                using FDTDacoustic<ValueType>::FDTDacoustic;
            };
        }
    }
}
