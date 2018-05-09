#include "FreeSurfaceElastic.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::~FreeSurfaceElastic(){};

/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model)
{

    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();

    auto temp = lama::eval<lama::DenseVector<ValueType>>( 2 * sWaveModulus - pWaveModulus );

    scaleHorizontalUpdate = 1.0 / pWaveModulus;
    scaleHorizontalUpdate *= temp;
    scaleHorizontalUpdate *= selectHorizontalUpdate;
}

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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Syy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative *= scaleHorizontalUpdate;

    Sxx += sumHorizonalDerivative;

    Syy *= setSurfaceZero;
}

/*! \brief Apply free surface condition during time stepping for 3D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 \param Szz Szz wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(scai::lama::Vector<ValueType> &sumHorizonalDerivative, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Syy, scai::lama::Vector<ValueType> &Szz)
{

    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative *= scaleHorizontalUpdate;

    Sxx += sumHorizonalDerivative;
    Szz += sumHorizonalDerivative;

    Syy *= setSurfaceZero;
}

/*! \brief Initialitation of the free surface
 *
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param NX Number of grid points in X-Direction
 \param NY Number of grid points in Y-Direction (Depth)
 \param NZ Number of grid points in Z-Direction
 \param DT Temporal Sampling
 \param DH Distance between grid points
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfPressure(NX, NY, NZ, dist);
    derivatives.calcDyfVelocity(NX, NY, NZ, dist);
    derivatives.calcDybPressure(NX, NY, NZ, dist);
    derivatives.calcDybVelocity(NX, NY, NZ, dist);
    derivatives.DybPressure *= DT / DH;
    derivatives.DybVelocity *= DT / DH;
    derivatives.DyfPressure *= DT / DH;
    derivatives.DyfVelocity *= DT / DH;
    derivatives.Dyb.purge();
    derivatives.Dyf.purge();

    selectHorizontalUpdate.setSameValue(dist, 0.0);
    setSurfaceZero.setSameValue(dist, 1.0);

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    /* Get write access to local part of scaleHorizontalUpdate */
    auto write_selectHorizontalUpdate = hostWriteAccess(selectHorizontalUpdate.getLocalValues());

    /* Get write access to local part of setSurfaceZero */
    auto write_setSurfaceZero = hostWriteAccess(setSurfaceZero.getLocalValues());

    KITGPI::Acquisition::Coordinates coordinateTransformation;

    IndexType rowGlobal;
    IndexType rowLocal;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        rowGlobal = read_localIndices[i];
        rowLocal = dist->global2local(rowGlobal);

        /* Determine if the current grid point is located on the surface */
        if (coordinateTransformation.locatedOnSurface(rowGlobal, NX, NY, NZ)) {

            /* Set horizontal update to 1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal] = 1.0;

            /* Set vector at surface to zero  */
            write_setSurfaceZero[rowLocal] = 0.0;
        }
    }

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<double>;
