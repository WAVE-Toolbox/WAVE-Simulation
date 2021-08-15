#pragma once

#include <scai/lama.hpp>

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
                FDTD3Dacoustic(){};
                //! Default destructor
                ~FDTD3Dacoustic(){};

            };
        }
    }
}
