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

/*! \brief Check wavefield
 *
 * Check wavefield for infs and nans
 *
 */
template <typename ValueType>
bool KITGPI::Wavefields::Wavefields<ValueType>::isFinite(scai::dmemo::DistributionPtr dist)
{
    bool result_isfinite=true;
    
    /* Get local owned global indices */
    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);
    
    auto read_VX = VX.getLocalValues();
    auto read_VZ = VZ.getLocalValues();

        IndexType localIndex=floor(ownedIndeces.size()/2);
        if (read_VX.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_VX[localIndex]),"Infinite or NaN value in VX wavefield!" << localIndex);
            if (isfinite(read_VX[localIndex])==false){
                result_isfinite=false;
            }
        }
        
        if (read_VZ.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_VZ[localIndex]),"Infinite or NaN value in VZ wavefield!" << localIndex);
            if (isfinite(read_VZ[localIndex])==false){
                result_isfinite=false;
            }
        }

    return(result_isfinite);
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

//! \brief Getter routine for vX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefVX()
{
    return (VX);
}

//! \brief Getter routine for vY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefVY()
{
    return (VY);
}

//! \brief Getter routine for vZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefVZ()
{
    return (VZ);
}

//! \brief Getter routine for Sxx wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSxx()
{
    return (Sxx);
}

//! \brief Getter routine for Syy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSyy()
{
    return (Syy);
}

//! \brief Getter routine for Szz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSzz()
{
    return (Szz);
}

//! \brief Getter routine for Syz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSyz()
{
    return (Syz);
}

//! \brief Getter routine for Sxz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSxz()
{
    return (Sxz);
}

//! \brief Getter routine for Sxy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefSxy()
{
    return (Sxy);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefP()
{
    return (P);
}

//! \brief Getter routine for Rxx Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRxx()
{
    return (Rxx);
}

//! \brief Getter routine for Ryy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRyy()
{
    return (Ryy);
}

//! \brief Getter routine for Rzz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRzz()
{
    return (Rzz);
}

//! \brief Getter routine for Ryz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRyz()
{
    return (Ryz);
}

//! \brief Getter routine for Rxz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRxz()
{
    return (Rxz);
}

//! \brief Getter routine for Rxy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRefRxy()
{
    return (Rxy);
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

template class KITGPI::Wavefields::Wavefields<float>;
template class KITGPI::Wavefields::Wavefields<double>;
