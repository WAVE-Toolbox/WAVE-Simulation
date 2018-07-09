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

            /*! \brief Class to create a Ricker sourcesignals
             *
             */
            template <typename ValueType>
            class Ricker : public SourceSignal<ValueType>
            {
              public:
                explicit Ricker(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~Ricker(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
