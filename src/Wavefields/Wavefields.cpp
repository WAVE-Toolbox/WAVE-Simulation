#include "Wavefields.hpp"
#include <math.h>

using namespace scai;

/*! \brief Reset a single wavefield to zero.
 \param vector Vector to be reset to 0
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::resetWavefield(scai::lama::DenseVector<ValueType> &vector)
{
    vector = 0;
}

/*! \brief Intitialisation of a single wavefield vector.
 *
 * This method will set the context, allocate the the wavefield and set the field to zero.
 *
 \param vector Vector to be set
 \param ctx Context pointer
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);

    resetWavefield(vector);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::Wavefields<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numWavefields)
{
    ValueType size = getMemoryWavefield(dist) / 1024 / 1024 * numWavefields;
    return size;
}

//! \brief calculate and return memory usage the of a Wavefield
/*!
 */
template <typename ValueType>
ValueType KITGPI::Wavefields::Wavefields<ValueType>::getMemoryWavefield(scai::dmemo::DistributionPtr dist)
{
    /* size of a wavefield is the size of a densevector = numGridpoints*size of Valuetype*/
    return (dist->getGlobalSize() * sizeof(ValueType));
}

/*! \brief Overloading = Operation
 *
 \param rhs Wavefield which is copied.
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::operator=(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Wavefield which is subtractet.
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::operator-=(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    minusAssign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Wavefield which is added.
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::operator+=(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::Wavefields<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::operator*=(ValueType rhs)
{
    timesAssign(rhs);
    return *this;
}

/*! \brief calculate a matrix to transform wavefield from XZY to YXZ
 * \param modelCoordinates coordinates of the original model
 * \param modelCoordinatesInversion coordinates of the averaged model
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::calcTransformMatrixXYZ(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    IndexType NX = modelCoordinates.getNX();
    IndexType NY = modelCoordinates.getNY();
    KITGPI::Acquisition::coordinate3D coordinateXZY;
    
    dmemo::DistributionPtr dist(transformMatrixXZY.getRowDistributionPtr());
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assemblyXZY;
    IndexType columnIndex;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        coordinateXZY = modelCoordinates.index2coordinate(ownedIndex);
        columnIndex = coordinateXZY.y + coordinateXZY.x * NY + coordinateXZY.z * NX * NY;
        assemblyXZY.push(ownedIndex, columnIndex, 1);
    }

    transformMatrixXZY.fillFromAssembly(assemblyXZY);
    transformMatrixYXZ.assignTranspose(transformMatrixXZY);
//     transformMatrixYXZ.writeToFile("wavefields/transformMatrixYXZ.mtx");
//     transformMatrixXZY.writeToFile("wavefields/transformMatrixXZY.mtx");
}

/*! \brief Initialize wavefield transform matrix to XYZ
 \param decomposeWavefieldType decomposeWavefieldType used to identify coordinate
 \param dist distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::initTransformMatrixXYZ(IndexType decomposeWavefieldType, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx)
{    
    if (decomposeWavefieldType != 0) {
        transformMatrixYXZ = lama::zero<SparseFormat>(dist, dist);
        transformMatrixYXZ.setContextPtr(ctx);
        transformMatrixXZY = lama::zero<SparseFormat>(dist, dist);
        transformMatrixXZY.setContextPtr(ctx);
    }
}

/*! \brief get transformMatrixYXZ
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::Wavefields<ValueType>::getTransformMatrixYXZ()
{
    return transformMatrixYXZ;
}

/*! \brief get transformMatrixXZY
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::Wavefields<ValueType>::getTransformMatrixXZY()
{
    return transformMatrixXZY;
}

template class KITGPI::Wavefields::Wavefields<float>;
template class KITGPI::Wavefields::Wavefields<double>;
