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

            //! \brief FDTD2Dtmem class
            template <typename ValueType>
            class FDTD3Demem : public SourceReceiverImplEM<ValueType>
            {
              public:
                //! Default constructor
                FDTD3Demem() = delete;
                //! Default destructor
                ~FDTD3Demem(){};

                using SourceReceiverImplEM<ValueType>::SourceReceiverImplEM;

              private:
                /* Temporary memory */
            };
        }
    }
}
