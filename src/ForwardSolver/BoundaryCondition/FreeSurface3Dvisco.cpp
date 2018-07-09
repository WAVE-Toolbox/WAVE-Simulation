#include "FreeSurface3Dvisco.hpp"

/*! \brief exchange the horizontal update at the free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonatlDerivative Sum of horizontal velocity updates
 \param vyy vertical velocity update
 \param Sxx Sxx wavefield (stress)
 \param Szz Szz wavefield (stress)
 \param Rxx Rxx wavefield (relaxation)
 \param Rzz Rzz wavefield (relaxation)
 \param DThalf Time discretisation DT/2
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<ValueType>::exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Szz, scai::lama::Vector<ValueType> &Rxx, scai::lama::Vector<ValueType> &Rzz, ValueType DThalf)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Undo Update of the relaxation parameter at the free surface */
    temp = selectHorizontalUpdate;
    temp *= Rxx;
    Sxx -= DThalf * temp;
    temp = selectHorizontalUpdate;
    temp *= Rzz;
    Szz -= DThalf * temp;

    /* Update the stress parameter at the free surface */
    /* vertical derivatives will be exchanged with horizontal derivatives (previous update will be undone) */
    temp = scaleStressHorizontalUpdate;
    temp *= sumHorizonatlDerivative;
    Sxx += temp;
    Szz += temp;

    temp = scaleStressVerticalUpdate;
    temp *= vyy;
    Sxx -= temp;
    Szz -= temp;

    /* Update relaxation parameter at the free surface */
    /* vertical derivatives will be exchanged with horizontal derivatives (previous update will be undone) */
    temp = scaleRelaxationHorizontalUpdate;
    temp *= sumHorizonatlDerivative;
    Rxx += temp;
    Rzz += temp;

    temp = scaleRelaxationVerticalUpdate;
    temp *= vyy;
    Rxx -= temp;
    Rzz -= temp;

    /* Apply update of the relaxation parameter to Sxx at the free surface */
    temp = selectHorizontalUpdate;
    temp *= Rxx;
    Sxx += DThalf * temp;
    temp = selectHorizontalUpdate;
    temp *= Rzz;
    Szz += DThalf * temp;
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<double>;
