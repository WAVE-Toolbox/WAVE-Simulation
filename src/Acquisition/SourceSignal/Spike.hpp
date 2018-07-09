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

            /*! \brief Class to create a Spike sourcesignals
             *
             */
            template <typename ValueType>
            class Spike : public SourceSignal<ValueType>
            {
              public:
                explicit Spike(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~Spike(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
