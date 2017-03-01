#include "FreeSurfaceAcoustic.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<ValueType>::~FreeSurfaceAcoustic(){};

/*! \brief Apply free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param p p wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<ValueType>::apply(scai::lama::Vector &p)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    /* Set the elements on the surface to zero */
    p.scale(setSurfaceZero);
}

/*! \brief Initialitation of the free surface
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfVelocity(NX, NY, NZ, dist);
    derivatives.DyfVelocity.scale(lama::Scalar(DT / DH));
    derivatives.Dyf.purge();

    /* Distributed vectors */
    setSurfaceZero.allocate(dist); // Vector to set elements on surface to zero
    setSurfaceZero = 1.0;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

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

            /* Set elements at the surface to zero */
            write_setSurfaceZero[rowLocal] = 0.0;
        }
    }
    read_localIndices.release();
    write_setSurfaceZero.release();

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<double>;