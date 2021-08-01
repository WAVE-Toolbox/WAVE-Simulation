#include "Hilbert.hpp"

using namespace scai;

/*! \brief set coefficient length
*
\param setCoefficientLength Length of coefficient vector
*/
template <typename ValueType>
void KITGPI::Hilbert::Hilbert<ValueType>::setCoefficientLength(IndexType const setCoefficientLength)
{
    kernelSize = setCoefficientLength;
}

template class KITGPI::Hilbert::Hilbert<double>;
template class KITGPI::Hilbert::Hilbert<float>;
