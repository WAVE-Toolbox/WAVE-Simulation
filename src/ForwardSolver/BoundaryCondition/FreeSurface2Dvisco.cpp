#include "FreeSurface2Dvisco.hpp"

/*! \brief Apply free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonatlDerivative Sum of horizontal velocity updates
 \param temp DenseVector which is used for temporary storage
 \param Sxx Sxx wavefield (stress)
 \param Syy Syy wavefield (stress)
 \param Rxx Rxx wavefield (relaxation)
 \param Ryy Ryy wavefield (relaxation)
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<ValueType>::apply(lama::Vector &sumHorizonatlDerivative, lama::Vector &temp, lama::Vector &Sxx, lama::Vector &Syy, lama::Vector &Rxx, lama::Vector &Ryy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Update the stress parameter at the free surface */
    temp = sumHorizonatlDerivative;
    temp.scale(scaleStressHorizontalUpdate);
    Sxx += temp; // Apply horizontal update

    Syy.scale(setSurfaceZero); // Set the free surface to zero

    /* Update relaxation parameter at the free surface */
    sumHorizonatlDerivative.scale(scaleRelaxationHorizontalUpdate);
    Rxx += sumHorizonatlDerivative; // Apply horizontal update

    Ryy.scale(setSurfaceZero); // Set the free surface to zero
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<double>;
