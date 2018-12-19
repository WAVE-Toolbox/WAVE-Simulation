#include "FreeSurface3Delastic.hpp"
using namespace scai;

/*! \brief Apply free surface condition during time stepping for 3D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param vyy vertical velocity update
 \param Scc Sxx wavefield
 \param Szz Szz wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Szz)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one */
    /* The previous Sxx update on the free surface will be undone */
    temp = scaleHorizontalUpdate;
    temp *= sumHorizonalDerivative;

    Sxx += temp;
    Szz += temp;

    temp = scaleVerticalUpdate;
    temp *= vyy;

    Sxx -= temp;
    Szz -= temp;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<double>;
