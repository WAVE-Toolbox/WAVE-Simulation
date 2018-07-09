#include "FreeSurface.hpp"

/*! \brief Getter method for active bool
 *
 *
 */
template <typename ValueType>
bool KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType>::getActive() const
{
    return (active);
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<double>;
