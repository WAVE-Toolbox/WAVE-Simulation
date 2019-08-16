#include "ABS3D.hpp"
using namespace scai;

/*! \brief Application of the damping boundary
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param v1 DenseVector to apply damping boundary
 \param v2 DenseVector to apply damping boundary
 \param v3 DenseVector to apply damping boundary
 \param v4 DenseVector to apply damping boundary
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::apply(lama::Vector<ValueType> &v1, lama::Vector<ValueType> &v2, lama::Vector<ValueType> &v3, lama::Vector<ValueType> &v4)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1 *= damping;
    v2 *= damping;
    v3 *= damping;
    v4 *= damping;
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
 \param v6 DenseVector to apply damping boundary
 \param v7 DenseVector to apply damping boundary
 \param v8 DenseVector to apply damping boundary
 \param v9 DenseVector to apply damping boundary
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::apply(
    lama::Vector<ValueType> &v1, lama::Vector<ValueType> &v2, lama::Vector<ValueType> &v3,
    lama::Vector<ValueType> &v4, lama::Vector<ValueType> &v5, lama::Vector<ValueType> &v6,
    lama::Vector<ValueType> &v7, lama::Vector<ValueType> &v8, lama::Vector<ValueType> &v9)
{

    SCAI_ASSERT_DEBUG(active, " ABS is not active ");

    v1 *= damping;
    v2 *= damping;
    v3 *= damping;
    v4 *= damping;
    v5 *= damping;
    v6 *= damping;
    v7 *= damping;
    v8 *= damping;
    v9 *= damping;
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
ValueType KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    /* Get local owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D coordinatedist;

    IndexType coordinateMin = 0;
    IndexType coordinatexzMin = 0;

    IndexType counter = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        coordinate = modelCoordinates.index2coordinate(ownedIndex);
        coordinatedist = modelCoordinates.edgeDistance(coordinate);

        coordinateMin = coordinatedist.min();

        if (useFreeSurface == 0) {
            if (coordinateMin < BoundaryWidth) {
                counter++;
            }

        } else {
            coordinatexzMin = !((coordinatedist.x) < (coordinatedist.z)) ? (coordinatedist.z) : (coordinatedist.x);
            if (coordinate.y < BoundaryWidth) {
                if ((coordinatedist.z < BoundaryWidth) || (coordinatedist.x < BoundaryWidth)) {
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
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the Damping Boundary...\n");

    active = true;

    /* Get local owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    /* Distributed vectors */
    damping.setSameValue(dist, 1.0); // sparse vector damping gets zero element 1.0

    // calculate damping function
    ValueType amp = 0;
    ValueType coeff[BoundaryWidth];
    ValueType a = 0;

    amp = 1.0 - DampingCoeff / 100.0;
    a = sqrt(-log(amp) / ((BoundaryWidth) * (BoundaryWidth)));

    for (IndexType j = 0; j < BoundaryWidth; j++) {
        coeff[j] = exp(-(a * a * (BoundaryWidth - j) * (BoundaryWidth - j)));
    }

    lama::VectorAssembly<ValueType> assembly;
    assembly.reserve(ownedIndeces.size());

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D coordinatedist;

    IndexType coordinateMin = 0;
    IndexType coordinatexzMin = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        coordinate = modelCoordinates.index2coordinate(ownedIndex);
        coordinatedist = modelCoordinates.edgeDistance(coordinate);

        coordinateMin = coordinatedist.min();

        if (useFreeSurface == 0) {
            if (coordinateMin < BoundaryWidth) {
                assembly.push(ownedIndex, coeff[coordinateMin]);
            }

        } else {
            coordinatexzMin = !((coordinatedist.x) < (coordinatedist.z)) ? (coordinatedist.z) : (coordinatedist.x);
            if (coordinate.y < BoundaryWidth) {
                if ((coordinatedist.z < BoundaryWidth) || (coordinatedist.x < BoundaryWidth)) {
                    assembly.push(ownedIndex, coeff[coordinatexzMin]);
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

template class KITGPI::ForwardSolver::BoundaryCondition::ABS3D<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::ABS3D<double>;
