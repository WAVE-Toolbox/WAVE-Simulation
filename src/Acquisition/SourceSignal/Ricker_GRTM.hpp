#pragma once

#include "SourceSignal.hpp"
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

namespace KITGPI
{

    namespace Acquisition
    {

        namespace SourceSignal
        {

            /*! \brief Class to create a Ricker_GRTM sourcesignals
             *
             */
            template <typename ValueType>
            class Ricker_GRTM : public SourceSignal<ValueType>
            {
              public:
                explicit Ricker_GRTM(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~Ricker_GRTM(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
