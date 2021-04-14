#include "FreeSurface2Dtmem.hpp"
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dtmem<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx)
{

    /* Apply horizontal update, which replaces the vertical one 
    * On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) (vxx) where mu = dielectricPermittivityEM and pi = velocivityEM
    * The original update,
    * sxx = pi * ( vxx+vyy) - 2mu *(vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx - 2mu / pi * vxx)
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx)  -  scaleVerticalUpdate*Hyy 
    * The update of szz is calculated the same way   
    */

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

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dtmem<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dtmem<double>;
