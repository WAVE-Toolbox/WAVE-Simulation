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

            /*! \brief Class to create a SinW sourcesignals
             *
             */
            template <typename ValueType>
            class SinW : public SourceSignal<ValueType>
            {
              public:
                explicit SinW(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);

                ///! Destructor
                ~SinW(){};

                void calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
            };
        }
    }
}
