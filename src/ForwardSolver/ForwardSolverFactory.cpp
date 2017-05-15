#include "ForwardSolverFactory.hpp"

template <typename ValueType>
typename KITGPI::ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr KITGPI::ForwardSolver::Factory<ValueType>::Create(std::string dimension, std::string type)
{

    // transform to lower cases
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0 || type.compare("sh") == 0, "Unkown type");

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("acoustic") == 0) {
        return ForwardSolverPtr(new FD2Dacoustic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("elastic") == 0) {
        return ForwardSolverPtr(new FD2Delastic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("visco") == 0) {
        return ForwardSolverPtr(new FD2Dvisco<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("sh") == 0) {
        return ForwardSolverPtr(new FD2Dsh<ValueType>);
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("acoustic") == 0) {
        return ForwardSolverPtr(new FD3Dacoustic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("elastic") == 0) {
        return ForwardSolverPtr(new FD3Delastic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("visco") == 0) {
        return ForwardSolverPtr(new FD3Dvisco<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return nullptr;
};

template class KITGPI::ForwardSolver::Factory<double>;
template class KITGPI::ForwardSolver::Factory<float>;
