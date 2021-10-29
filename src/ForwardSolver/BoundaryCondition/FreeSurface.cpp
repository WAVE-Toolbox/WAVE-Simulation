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

template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType>::setSurfaceZero(scai::lama::Vector<ValueType> &vector)
{
    /* this function is actually unnecessary because syy will be set implicitly to zero during the update of the velocities. Nevertheless for further use of the wavefields (i.e. in FWI) its better to
    set the values to zero because they have arbitrary< high values which disturb the gradient */

    vector *= setZeroFreeSurface;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<double>;
