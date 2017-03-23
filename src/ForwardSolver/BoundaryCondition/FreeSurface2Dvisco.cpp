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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<ValueType>::apply(scai::lama::Vector &sumHorizonatlDerivative, scai::lama::Vector &temp, scai::lama::Vector &Sxx, scai::lama::Vector &Syy, scai::lama::Vector &Rxx, scai::lama::Vector &Ryy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Update the stress parameter at the free surface */
    temp = sumHorizonatlDerivative;
    temp *= scaleStressHorizontalUpdate;
    Sxx += temp; // Apply horizontal update

    Syy *= setSurfaceZero; // Set the free surface to zero

    /* Update relaxation parameter at the free surface */
    sumHorizonatlDerivative *= scaleRelaxationHorizontalUpdate;
    Rxx += sumHorizonatlDerivative; // Apply horizontal update

    Ryy *= setSurfaceZero; // Set the free surface to zero
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco<double>;
