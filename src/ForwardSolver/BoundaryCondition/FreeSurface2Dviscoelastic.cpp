#include "FreeSurface2Dviscoelastic.hpp"

/*! \brief exchange the horizontal update at the free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonatlDerivative Sum of horizontal velocity updates
 \param vyy vertical velocity update
 \param Sxx Sxx wavefield (stress)
 \param Rxx Rxx wavefield (relaxation)
 \param DThalf Time discretisation DT/2
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dviscoelastic<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, std::vector<scai::lama::DenseVector<ValueType>> &Rxx, ValueType DThalf)
{
    /* On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) * (vxx+vzz) where mu = sWaveModulus * (1+L*tauS) and pi = pWaveModulus * (1+L*tauP)
    * The original update,
    * sxx = pi * ( vxx+vyy+vzz ) - 2mu *(vzz+vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx + vzz - 2mu / pi * (vxx+vzz))
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx+vzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx+vzz)  - scaleVerticalUpdate*Vyy 
    * the szz calculation is done analogous */

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Undo Update of the relaxation parameter at the free surface */
    for (int l=0; l<numRelaxationMechanisms; l++) {
        temp = selectFreeSurface;
        temp *= Rxx[l];
        Sxx -= DThalf * temp;
    }

    /* Update the stress parameter at the free surface */
    /* vertical derivatives will be exchanged with horizontal derivatives (previous update will be undone) */
    temp = scaleStressHorizontalUpdate;
    temp *= sumHorizonatlDerivative;
    Sxx += temp; // Apply horizontal update

    temp = scaleStressVerticalUpdate;
    temp *= vyy;
    Sxx -= temp;

    /* Update relaxation parameter at the free surface */
    /* vertical derivatives will be exchanged with horizontal derivatives (previous update will be undone) */
    for (int l=0; l<numRelaxationMechanisms; l++) {
        temp = scaleRelaxationHorizontalUpdate[l];
        temp *= sumHorizonatlDerivative;
        Rxx[l] += temp; // Apply horizontal update

        temp = scaleRelaxationVerticalUpdate[l];
        temp *= vyy;
        Rxx[l] -= temp;

        /* Apply update of the relaxation parameter to Sxx at the free surface */
        temp = selectFreeSurface;
        temp *= Rxx[l];
        Sxx += DThalf * temp;
    }
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dviscoelastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dviscoelastic<double>;
