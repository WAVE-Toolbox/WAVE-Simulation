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
            class FDTD2Demem : public SourceReceiverImplEM<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Demem() = delete;
                //! Default destructor
                ~FDTD2Demem(){};

                using SourceReceiverImplEM<ValueType>::SourceReceiverImplEM;

                void applySourceEZ(Acquisition::SeismogramEM<ValueType> const & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EZ sources can not be implemented in EMEM modeling")};
                void gatherSeismogramEZ(Acquisition::SeismogramEM<ValueType> & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EZ receivers can not be implementedin EMEM modeling")};

              private:
                /* Temporary memory */
            };
        }
    }
}
