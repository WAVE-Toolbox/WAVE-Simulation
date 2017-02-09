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

    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &sWaveModulus = model.getSWaveModulus();

    lama::DenseVector<ValueType> temp(sWaveModulus.getDistributionPtr());

    temp = 2 * sWaveModulus - pWaveModulus;

    scaleHorizontalUpdate = pWaveModulus;
    scaleHorizontalUpdate.invert();
    scaleHorizontalUpdate.scale(temp);
    scaleHorizontalUpdate.scale(selectHorizontalUpdate);
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(scai::lama::Vector &sumHorizonalDerivative, scai::lama::Vector &Sxx, scai::lama::Vector &Syy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);

    Sxx += sumHorizonalDerivative;

    Syy.scale(setSurfaceZero);
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(scai::lama::Vector &sumHorizonalDerivative, scai::lama::Vector &Sxx, scai::lama::Vector &Syy, scai::lama::Vector &Szz)
{

    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);

    Sxx += sumHorizonalDerivative;
    Szz += sumHorizonalDerivative;

    Syy.scale(setSurfaceZero);
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
    derivatives.DybPressure.scale(lama::Scalar(DT / DH));
    derivatives.DybVelocity.scale(lama::Scalar(DT / DH));
    derivatives.DyfPressure.scale(lama::Scalar(DT / DH));
    derivatives.DyfVelocity.scale(lama::Scalar(DT / DH));
    derivatives.Dyb.purge();
    derivatives.Dyf.purge();

    selectHorizontalUpdate.allocate(dist);
    selectHorizontalUpdate = 0.0;

    setSurfaceZero.allocate(dist);
    setSurfaceZero = 1.0;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    /* Get write access to local part of scaleHorizontalUpdate */
    utilskernel::LArray<ValueType> *selectHorizontalUpdate_LA = &selectHorizontalUpdate.getLocalValues();
    hmemo::WriteAccess<ValueType> write_selectHorizontalUpdate(*selectHorizontalUpdate_LA);

    /* Get write access to local part of setSurfaceZero */
    utilskernel::LArray<ValueType> *setSurfaceZero_LA = &setSurfaceZero.getLocalValues();
    hmemo::WriteAccess<ValueType> write_setSurfaceZero(*setSurfaceZero_LA);

    KITGPI::Acquisition::Coordinates<ValueType> coordinateTransformation;

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
    read_localIndices.release();
    write_selectHorizontalUpdate.release();
    write_setSurfaceZero.release();
    HOST_PRINT(dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<double>;
