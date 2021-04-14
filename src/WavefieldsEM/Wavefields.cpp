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

template class KITGPI::Wavefields::WavefieldsEM<float>;
template class KITGPI::Wavefields::WavefieldsEM<double>;
