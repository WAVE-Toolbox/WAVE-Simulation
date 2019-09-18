#include "FreeSurface2Delastic.hpp"
using namespace scai;

/*! \brief exchange Sxx at the free surface during time stepping for 2D elaastic simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param vyy vyy velocity updates
 \param Sxx Sxx wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one */
    /* The previous Sxx update on the free surface will be undone */

    // scaleHorizontalUpdate is a sparse vector with non zeors at the free surface
    temp = scaleHorizontalUpdate;
    temp *= sumHorizonalDerivative;

    Sxx += temp;
    
    // scaleVerticalUpdate is a sparse vector with non zeors at the free surface
    temp = scaleVerticalUpdate;
    temp *= vyy;

    Sxx -= temp;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<double>;
