#include "Wavefields.hpp"
#include <math.h>

using namespace scai;

/*! \brief Reset a single wavefield to zero.
 \param vector Vector to be reset to 0
 */
template <typename ValueType>
void KITGPI::Wavefields::WavefieldsEM<ValueType>::resetWavefield(scai::lama::DenseVector<ValueType> &vector)
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
void KITGPI::Wavefields::WavefieldsEM<ValueType>::initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);

    resetWavefield(vector);
}

/*! \brief Check wavefield
 *
 * Check wavefield for infs and nans
 *
 */
template <typename ValueType>
bool KITGPI::Wavefields::WavefieldsEM<ValueType>::isFinite(scai::dmemo::DistributionPtr dist)
{
    bool result_isfinite=true;
    
    /* Get local owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);
    
    auto read_EX = EX.getLocalValues();
    auto read_EY = EY.getLocalValues();
    auto read_EZ = EZ.getLocalValues();

    for (IndexType localIndex=0;localIndex<ownedIndeces.size();localIndex++) {
        if (read_EX.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_EX[localIndex]),"Infinite or NaN value in EX wavefield!" << localIndex);
            if (isfinite(read_EX[localIndex])==false){
                result_isfinite=false;
            }
        }
        
        if (read_EY.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_EY[localIndex]),"Infinite or NaN value in EY wavefield!" << localIndex);
            if (isfinite(read_EY[localIndex])==false){
                result_isfinite=false;
            }
        }
        
        if (read_EZ.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_EZ[localIndex]),"Infinite or NaN value in EZ wavefield!" << localIndex);
            if (isfinite(read_EZ[localIndex])==false){
                result_isfinite=false;
            }
        }
    }

    return(result_isfinite);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::WavefieldsEM<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numWavefields)
{
    ValueType size = getMemoryWavefield(dist) / 1024 / 1024 * numWavefields;
    return size;
}

//! \brief calculate and return memory usage the of a Wavefield
/*!
 */
template <typename ValueType>
ValueType KITGPI::Wavefields::WavefieldsEM<ValueType>::getMemoryWavefield(scai::dmemo::DistributionPtr dist)
{
    /* size of a wavefield is the size of a densevector = numGridpoints*size of Valuetype*/
    return (dist->getGlobalSize() * sizeof(ValueType));
}

//! \brief Getter routine for HX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefHX()
{
    return (HX);
}

//! \brief Getter routine for HY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefHY()
{
    return (HY);
}

//! \brief Getter routine for HZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefHZ()
{
    return (HZ);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEX()
{
    return (EX);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEY()
{
    return (EY);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEZ()
{
    return (EZ);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEXup()
{
    return (EXup);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEYup()
{
    return (EYup);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEZup()
{
    return (EZup);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEXdown()
{
    return (EXdown);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEYdown()
{
    return (EYdown);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEZdown()
{
    return (EZdown);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEXleft()
{
    return (EXleft);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEYleft()
{
    return (EYleft);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEZleft()
{
    return (EZleft);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEXright()
{
    return (EXright);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEYright()
{
    return (EYright);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefEZright()
{
    return (EZright);
}

//! \brief Getter routine for RX Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRX()
{
    return (RX);
}

//! \brief Getter routine for RY Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRY()
{
    return (RY);
}

//! \brief Getter routine for RZ Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRZ()
{
    return (RZ);
}

/*! \brief Overloading = Operation
 *
 \param rhs Wavefield which is copied.
 */
template <typename ValueType>
KITGPI::Wavefields::WavefieldsEM<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::operator=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Wavefield which is subtractet.
 */
template <typename ValueType>
KITGPI::Wavefields::WavefieldsEM<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::operator-=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    minusAssign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Wavefield which is added.
 */
template <typename ValueType>
KITGPI::Wavefields::WavefieldsEM<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::operator+=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::WavefieldsEM<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::operator*=(ValueType rhs)
{
    timesAssign(rhs);
    return *this;
}

/*! \brief calculate a matrix to transform wavefield from XZY to YXZ
 * \param modelCoordinates coordinates of the original model
 * \param modelCoordinatesInversion coordinates of the averaged model
 */
template <typename ValueType>
void KITGPI::Wavefields::WavefieldsEM<ValueType>::calcTransformMatrixXYZ(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates)
{
    // this function is only valid for grid partitioning with useVariableGrid=0
    IndexType NX = modelCoordinates.getNX();
    IndexType NY = modelCoordinates.getNY();
    IndexType NZ = modelCoordinates.getNZ();
    KITGPI::Acquisition::coordinate3D coordinateXZY;
    
    dmemo::DistributionPtr dist(transformMatrixXZY.getRowDistributionPtr());
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assemblyXZY;
    lama::MatrixAssembly<ValueType> assemblyX;
    lama::MatrixAssembly<ValueType> assemblyY;
    IndexType columnIndex;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        coordinateXZY = modelCoordinates.index2coordinate(ownedIndex);
        columnIndex = coordinateXZY.y + coordinateXZY.x * NY + coordinateXZY.z * NX * NY;
        assemblyXZY.push(ownedIndex, columnIndex, 1);
        if (coordinateXZY.x % 2 == 0) {
            if (coordinateXZY.x % NX == 0) {
                assemblyX.push(ownedIndex, ownedIndex+1, 1); 
            } else if (coordinateXZY.x % NX == NX-1) {
                assemblyX.push(ownedIndex, ownedIndex-1, 1); 
            } else {
                assemblyX.push(ownedIndex, ownedIndex+1, 0.5); 
                assemblyX.push(ownedIndex, ownedIndex-1, 0.5);
            }
        } else {
            assemblyX.push(ownedIndex, ownedIndex, 1);
        }
        if (coordinateXZY.y % 2 == 0) {
            if (coordinateXZY.y % NY == 0) {
                assemblyY.push(ownedIndex, ownedIndex+NX*NZ, 1); 
            } else if (coordinateXZY.y % NY == NX-1) {
                assemblyY.push(ownedIndex, ownedIndex-NX*NZ, 1); 
            } else {
                assemblyY.push(ownedIndex, ownedIndex+NX*NZ, 0.5); 
                assemblyY.push(ownedIndex, ownedIndex-NX*NZ, 0.5);
            }
        } else {
            assemblyY.push(ownedIndex, ownedIndex, 1);
        }
    }

    transformMatrixXZY.fillFromAssembly(assemblyXZY);
    transformMatrixYXZ.assignTranspose(transformMatrixXZY);
    transformMatrixX.fillFromAssembly(assemblyX);
    transformMatrixY.fillFromAssembly(assemblyY);
//     transformMatrixYXZ.writeToFile("wavefields/transformMatrixYXZ.mtx");
//     transformMatrixXZY.writeToFile("wavefields/transformMatrixXZY.mtx");
//     transformMatrixX.writeToFile("wavefields/transformMatrixX.mtx");
}

/*! \brief Initialize wavefield transform matrix to XYZ
 \param decomposeType decomposeType used to identify coordinate
 \param dist distribution
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Wavefields::WavefieldsEM<ValueType>::initTransformMatrixXYZ(IndexType decomposeType, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx)
{    
    if (decomposeType != 0) {
        transformMatrixYXZ = lama::zero<SparseFormat>(dist, dist);
        transformMatrixYXZ.setContextPtr(ctx);
        transformMatrixXZY = lama::zero<SparseFormat>(dist, dist);
        transformMatrixXZY.setContextPtr(ctx);
        transformMatrixX = lama::zero<SparseFormat>(dist, dist);
        transformMatrixX.setContextPtr(ctx);
        transformMatrixY = lama::zero<SparseFormat>(dist, dist);
        transformMatrixY.setContextPtr(ctx);
    }
}

/*! \brief get transformMatrixYXZ
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::WavefieldsEM<ValueType>::getTransformMatrixYXZ()
{
    return transformMatrixYXZ;
}

/*! \brief get transformMatrixXZY
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::WavefieldsEM<ValueType>::getTransformMatrixXZY()
{
    return transformMatrixXZY;
}

/*! \brief get transformMatrixX
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::WavefieldsEM<ValueType>::getTransformMatrixX()
{
    return transformMatrixX;
}
/*! \brief get transformMatrixY
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> const &KITGPI::Wavefields::WavefieldsEM<ValueType>::getTransformMatrixY()
{
    return transformMatrixY;
}


template class KITGPI::Wavefields::WavefieldsEM<float>;
template class KITGPI::Wavefields::WavefieldsEM<double>;
