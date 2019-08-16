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

//! \brief estimate memory for the absorbing boundary frame
/*!
 *
 \param BoundaryWidth Width of damping boundary
 \param useFreeSurface Indicator which free surface is in use
  \param dist Distribution of the wavefield
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    /* Get local owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    IndexType counter = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        Acquisition::coordinate3D coordinatedist = modelCoordinates.edgeDistance(coordinate);

        IndexType temp = 0;
        if (coordinatedist.x < coordinatedist.y) {
            temp = coordinatedist.x;
        } else {
            temp = coordinatedist.y;
        }
        IndexType coordinateMin = temp;
        if (useFreeSurface == 0) {
            if (coordinateMin < BoundaryWidth) {
                counter++;
            }

        } else {
            if (coordinate.y < BoundaryWidth) {
                if ((coordinatedist.x < BoundaryWidth)) {
                    counter++;
                }
            } else if (coordinateMin < BoundaryWidth) {
                counter++;
            }
        }
    }

    IndexType numPartitions = dist->getNumPartitions();
    ValueType mega = 1024 * 1024;
    ValueType size = counter * sizeof(ValueType) / mega;
    HOST_PRINT(dist->getCommunicatorPtr(), " -  Absorbing Frame Vector  \t" << size << " / " << size / numPartitions << " MB\n");
    return size;
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
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the Damping Boundary...\n");

    active = true;

    /* Get owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);
    /* Distributed vectors */
    /*allocate Sparse Vector and temporary Dense Vector*/
    damping.setSameValue(dist, 1.0); // sparse vector damping gets zero element 1.0

    lama::VectorAssembly<ValueType> assembly;
    assembly.reserve(ownedIndeces.size());

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

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        coordinate = modelCoordinates.index2coordinate(ownedIndex);
        coordinatedist = modelCoordinates.edgeDistance(coordinate);

        temp = 0;
        if (coordinatedist.x < coordinatedist.y) {
            temp = coordinatedist.x;
        } else {
            temp = coordinatedist.y;
        }
        coordinateMin = temp;
        if (useFreeSurface == 0) {
            if (coordinateMin < BoundaryWidth) {
                assembly.push(ownedIndex, coeff[coordinateMin]);
            }

        } else {
            if (coordinate.y < BoundaryWidth) {
                if ((coordinatedist.x < BoundaryWidth)) {
                    assembly.push(ownedIndex, coeff[coordinatedist.x]);
                }
            } else if (coordinateMin < BoundaryWidth) {
                assembly.push(ownedIndex, coeff[coordinateMin]);
            }
        }
    }

    damping.setContextPtr(ctx);
    damping.fillFromAssembly(assembly);

    HOST_PRINT(comm, "Finished with initialization of the Damping Boundary!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::ABS2D<float>;
