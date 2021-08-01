#include "FreeSurface3Demem.hpp"
using namespace scai;

/*! \brief Apply free surface condition during time stepping for 3D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param vyy vertical velocity update
 \param Sxx Sxx wavefield
 \param Szz Szz wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Demem<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Szz)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one 
    * On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) (vxx+hzz) where mu = dielectricPermittivity and pi = velocivityEM
    * The original update,
    * sxx = pi * ( vxx+vyy+hzz ) - 2mu *(hzz+vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx + hzz - 2mu / pi * (vxx+hzz))
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx + hzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx + hzz)  -  scaleVerticalUpdate*Hyy 
    * The update of szz is calculated the same way   
    */
    
    
    // scaleHorizontalUpdate is a sparse vector with non zeors at the free surface
    temp = scaleHorizontalUpdate;
    temp *= sumHorizonalDerivative;

    Sxx += temp;
    Szz += temp;

    // scaleVerticalUpdate is a sparse vector with non zeors at the free surface
    temp = scaleVerticalUpdate;
    temp *= vyy;

    Sxx -= temp;
    Szz -= temp;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Demem<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Demem<double>;
