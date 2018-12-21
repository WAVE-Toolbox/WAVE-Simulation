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
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1 *= damping;
    v2 *= damping;
    v3 *= damping;
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
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4, scai::lama::Vector<ValueType> &v5)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1 *= damping;
    v2 *= damping;
    v3 *= damping;
    v4 *= damping;
    v5 *= damping;
}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param BoundaryWidth Width of damping boundary
 \param DampingCoeff Damping coefficient
 \param useFreeSurface Indicator which free surface is in use
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx,Acquisition::Coordinates const &modelCoordinates, IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the Damping Boundary...\n");

    active = true;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);

    IndexType numLocalIndices = localIndices.size(); // Number of local indices

    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp;                             // Temporary storage, so we do not have to access the array

    /* Distributed vectors */
    /*allocate Sparse Vector and temporary Dense Vector*/
    damping.setSameValue(dist, 1.0); // sparse vector damping gets zero element 1.0
    lama::DenseVector<ValueType> damping_temp(dist, 1.0);

    /* Get write access to local part of damping_temp */
    hmemo::HArray<ValueType> *damping_LA = &damping_temp.getLocalValues();
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

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D coordinatedist;

    IndexType coordinateMin = 0;
    IndexType temp;

    /* Set the values into the indice arrays and the value array */
    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = modelCoordinates.index2coordinate(read_localIndices_temp);
        coordinatedist = modelCoordinates.edgeDistance(coordinate);

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

        if (useFreeSurface > 0) {
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
    damping = damping_temp;
    HOST_PRINT(comm, "Finished with initialization of the Damping Boundary!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<float>;
