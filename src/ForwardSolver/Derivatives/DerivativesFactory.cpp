#include "DerivativesFactory.hpp"
using namespace scai;

template <typename ValueType>
typename KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr KITGPI::ForwardSolver::Derivatives::Factory<ValueType>::Create(std::string dimension)
{

    // transform to lower cases
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("3d") == 0 || dimension.compare("2d") == 0, "Unkown dimension");

    if (dimension.compare("2d") == 0) {
        return DerivativesPtr(new FDTD2D<ValueType>);
    }
    if (dimension.compare("3d") == 0) {
        return DerivativesPtr(new FDTD3D<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return DerivativesPtr();
}

template class KITGPI::ForwardSolver::Derivatives::Factory<float>;
template class KITGPI::ForwardSolver::Derivatives::Factory<double>;
