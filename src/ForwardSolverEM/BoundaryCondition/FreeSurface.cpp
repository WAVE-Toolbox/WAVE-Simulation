#include "FreeSurface.hpp"

/*! \brief Getter method for active bool
 *
 *
 */
template <typename ValueType>
bool KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<ValueType>::getActive() const
{
    return (active);
}

template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<ValueType>::setSurfaceZero(scai::lama::Vector<ValueType> &vector)
{
    /* this function is actually unnecessary because syy will be set implicitly to zero during the update of the velocities. Nethertheless for further use of the wavefields (i.e. in FWI) its better to
    set the values to zero because they have abitrary< high values which disturb the gradient */

    vector *= setZeroFreeSurface;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<double>;
