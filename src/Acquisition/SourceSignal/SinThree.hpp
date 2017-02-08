#pragma once

#include "SourceSignal.hpp"
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

using namespace scai;

namespace KITGPI
{

    namespace Acquisition
    {

        namespace SourceSignal
        {

            /*! \brief Class to create a SinThree sourcesignals
             *
             */
            template <typename ValueType>
            class SinThree : public SourceSignal<ValueType>
            {
              public:
                explicit SinThree(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~SinThree(){};

                void calc(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
