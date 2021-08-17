#pragma once

#include <scai/lama.hpp>

#include "SourceReceiverImplEM.hpp"

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
                
                void applySourceEZ(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EZ sources can not be implemented in EMEM modeling")};
                void gatherSeismogramEZ(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::Wavefields<ValueType> & /*wavefieldIN*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EZ receivers can not be implemented in EMEM modeling")};

              private:
                /* Temporary memory */
            };
        }
    }
}
