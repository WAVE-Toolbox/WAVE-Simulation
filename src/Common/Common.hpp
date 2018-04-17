#pragma once

#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>

namespace KITGPI
{
    //! \brief Common namespace
    namespace Common
    {
        template<typename ValueType, typename IndexType> void searchAndReplace(scai::lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, IndexType compareType);
    }
}

# include "./Common.tpp"