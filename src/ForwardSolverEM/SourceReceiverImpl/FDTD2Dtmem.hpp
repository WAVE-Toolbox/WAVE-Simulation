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
            class FDTD2Dtmem : public SourceReceiverImplEM<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Dtmem() = delete;
                //! Default destructor
                ~FDTD2Dtmem(){};

                using SourceReceiverImplEM<ValueType>::SourceReceiverImplEM;

                void applySourceEX(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EX sources can not be implemented in 2D tmem modeling")};
                void gatherSeismogramEX(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EX receivers can not be implemented in 2D tmem modeling")};
                void applySourceEY(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EY sources can not be implemented in 2D tmem modeling")};
                void gatherSeismogramEY(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("EY receivers can not be implemented in 2D tmem modeling")};
                void applySourceHZ(Acquisition::Seismogram<ValueType> const & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("HZ sources can not be implemented in 2D tmem modeling")};
                void gatherSeismogramHZ(Acquisition::Seismogram<ValueType> & /*seismo*/, Wavefields::WavefieldsEM<ValueType> & /*wavefield*/, scai::IndexType /*t*/) override{COMMON_THROWEXCEPTION("HZ receivers can not be implemented in 2D tmem modeling")};

              private:
                /* Temporary memory */
            };
        }
    }
}
