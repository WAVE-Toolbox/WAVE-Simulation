
#pragma once

#include <scai/lama.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

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
            class FDTD2Delastic : public SourceReceiverImpl<ValueType>
            {
              public:
                //! Default constructor
                FDTD2Delastic() = delete;
                //! Default destructor
                ~FDTD2Delastic(){};

                using SourceReceiverImpl<ValueType>::SourceReceiverImpl;

                void applySourcePressure(Acquisition::Seismogram<ValueType> const &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t) override;
                void gatherSeismogramPressure(Acquisition::Seismogram<ValueType> &seismo, Wavefields::Wavefields<ValueType> &wavefield, IndexType t) override;

              private:
                /* Temporary memory */
                using SourceReceiverImpl<ValueType>::applySource_samplesPressure;
                using SourceReceiverImpl<ValueType>::gatherSeismogram_samplesPressure;
            };
        }
    }
}
