#include "FreeSurface2Delastic.hpp"
using namespace scai;

/*! \brief Apply free surface condition during time stepping for 2D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<ValueType>::apply(scai::lama::Vector &sumHorizonalDerivative, scai::lama::Vector &Sxx, scai::lama::Vector &Syy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative *= scaleHorizontalUpdate;

    Sxx += sumHorizonalDerivative;

    Syy *= setSurfaceZero;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<double>;
