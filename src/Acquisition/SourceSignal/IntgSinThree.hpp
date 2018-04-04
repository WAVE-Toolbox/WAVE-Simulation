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

            /*! \brief Class to create a IntgSinThree sourcesignals
             *
             */
            template <typename ValueType>
            class IntgSinThree : public SourceSignal<ValueType>
            {
              public:
                explicit IntgSinThree(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~IntgSinThree(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
