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

            /*! \brief Class to create a Berlage sourcesignals
             *
             */
            template <typename ValueType>
            class Berlage : public SourceSignal<ValueType>
            {
              public:
                explicit Berlage(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~Berlage(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, scai::IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
