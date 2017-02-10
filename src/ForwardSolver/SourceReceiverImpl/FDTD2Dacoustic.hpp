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

                using FDTDacoustic<ValueType>::FDTDacoustic;
            };
        }
    }
}
