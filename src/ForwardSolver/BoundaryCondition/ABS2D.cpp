#include "ABS2D.hpp"
using namespace scai;

/*! \brief Application of the damping boundary
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param v1 DenseVector to apply damping boundary
 \param v2 DenseVector to apply damping boundary
 \param v3 DenseVector to apply damping boundary
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
}

/*! \brief Application of the damping boundary
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param v1 DenseVector to apply damping boundary
 \param v2 DenseVector to apply damping boundary
 \param v3 DenseVector to apply damping boundary
 \param v4 DenseVector to apply damping boundary
 \param v5 DenseVector to apply damping boundary
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3, scai::lama::Vector &v4, scai::lama::Vector &v5)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
    v4.scale(damping);
    v5.scale(damping);
}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param BoundaryWidth Width of damping boundary
 \param DampingCoeff Damping coefficient
 \param useFreeSurface Bool if free surface is in use
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the Damping Boundary...\n");

    active = true;

    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);

    IndexType numLocalIndices = localIndices.size(); // Number of local indices

    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp;                             // Temporary storage, so we do not have to access the array

    /* Distributed vectors */
    damping.allocate(dist); // Vector to set elements on surface to zero
    damping = 1.0;

    /* Get write access to local part of setSurfaceZero */
    utilskernel::LArray<ValueType> *damping_LA = &damping.getLocalValues();
    hmemo::WriteAccess<ValueType> write_damping(*damping_LA);

    // calculate damping function
    ValueType amp = 0;
    ValueType coeff[BoundaryWidth];
    ValueType a = 0;

    amp = 1.0 - DampingCoeff / 100.0;
    a = sqrt(-log(amp) / ((BoundaryWidth) * (BoundaryWidth)));

    for (IndexType j = 0; j < BoundaryWidth; j++) {
        coeff[j] = exp(-(a * a * (BoundaryWidth - j) * (BoundaryWidth - j)));
    }

    Acquisition::Coordinates<ValueType> coordTransform;
    SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(2, 100, 100, 100).x == 2, "")
    SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(102, 100, 100, 100).y == 1, "")
    SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(2, 100, 100, 1).z == 0, "")

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D coordinatedist;

    IndexType coordinateMin = 0;
    IndexType temp;

    /* Set the values into the indice arrays and the value array */
    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ);
        coordinatedist = coordTransform.edgeDistance(coordinate, NX, NY, NZ);

        temp = 0;
        if (coordinatedist.x < coordinatedist.y) {
            temp = coordinatedist.x;
        } else {
            temp = coordinatedist.y;
        }
        coordinateMin = temp;
        if (coordinateMin < BoundaryWidth) {
            write_damping[i] = coeff[coordinateMin];
        }

        if (useFreeSurface) {
            if (coordinate.y < BoundaryWidth) {
                write_damping[i] = 1.0;

                if ((coordinatedist.x < BoundaryWidth)) {
                    write_damping[i] = coeff[coordinatedist.x];
                }
            }
        }
    }

    /* Release all read and write access */
    read_localIndices.release();
    write_damping.release();

    damping.setContextPtr(ctx);

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished with initialization of the Damping Boundary!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<float>;
