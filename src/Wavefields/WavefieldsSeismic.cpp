#include "WavefieldsSeismic.hpp"
#include <math.h>

using namespace scai;

/*! \brief Check wavefield
 *
 * Check wavefield for infs and nans
 *
 */
template <typename ValueType>
bool KITGPI::Wavefields::WavefieldsSeismic<ValueType>::isFinite(scai::dmemo::DistributionPtr dist)
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

//! \brief Getter routine for vX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVX()
{
    return (VX);
}

//! \brief Getter routine for vY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVY()
{
    return (VY);
}

//! \brief Getter routine for vZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVZ()
{
    return (VZ);
}

//! \brief Getter routine for Sxx wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSxx()
{
    return (Sxx);
}

//! \brief Getter routine for Syy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSyy()
{
    return (Syy);
}

//! \brief Getter routine for Szz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSzz()
{
    return (Szz);
}

//! \brief Getter routine for Syz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSyz()
{
    return (Syz);
}

//! \brief Getter routine for Sxz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSxz()
{
    return (Sxz);
}

//! \brief Getter routine for Sxy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefSxy()
{
    return (Sxy);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefP()
{
    return (P);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefPup()
{
    return (Pup);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefPdown()
{
    return (Pdown);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefPleft()
{
    return (Pleft);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefPright()
{
    return (Pright);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVXup()
{
    return (VXup);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVXdown()
{
    return (VXdown);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVXleft()
{
    return (VXleft);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVXright()
{
    return (VXright);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVYup()
{
    return (VYup);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVYdown()
{
    return (VYdown);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVYleft()
{
    return (VYleft);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVYright()
{
    return (VYright);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVZup()
{
    return (VZup);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVZdown()
{
    return (VZdown);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVZleft()
{
    return (VZleft);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefVZright()
{
    return (VZright);
}

//! \brief Getter routine for Rxx Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRxx()
{
    return (Rxx);
}

//! \brief Getter routine for Ryy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRyy()
{
    return (Ryy);
}

//! \brief Getter routine for Rzz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRzz()
{
    return (Rzz);
}

//! \brief Getter routine for Ryz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRyz()
{
    return (Ryz);
}

//! \brief Getter routine for Rxz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRxz()
{
    return (Rxz);
}

//! \brief Getter routine for Rxy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRxy()
{
    return (Rxy);
}



/* EM */
//! \brief Getter routine for HX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefHX()
{
    COMMON_THROWEXCEPTION("There is no HX in an Seismic modelling")
    return (HX);
}

//! \brief Getter routine for HY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefHY()
{
    COMMON_THROWEXCEPTION("There is no HY in an Seismic modelling")
    return (HY);
}

//! \brief Getter routine for HZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefHZ()
{
    COMMON_THROWEXCEPTION("There is no HZ in an Seismic modelling")
    return (HZ);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEX()
{
    COMMON_THROWEXCEPTION("There is no EX in an Seismic modelling")
    return (EX);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEY()
{
    COMMON_THROWEXCEPTION("There is no EY in an Seismic modelling")
    return (EY);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEZ()
{
    COMMON_THROWEXCEPTION("There is no EZ in an Seismic modelling")
    return (EZ);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEXup()
{
    COMMON_THROWEXCEPTION("There is no EXup in an Seismic modelling")
    return (EXup);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEYup()
{
    COMMON_THROWEXCEPTION("There is no EYup in an Seismic modelling")
    return (EYup);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEZup()
{
    COMMON_THROWEXCEPTION("There is no EZup in an Seismic modelling")
    return (EZup);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEXdown()
{
    COMMON_THROWEXCEPTION("There is no EXdown in an Seismic modelling")
    return (EXdown);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEYdown()
{
    COMMON_THROWEXCEPTION("There is no EYdown in an Seismic modelling")
    return (EYdown);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEZdown()
{
    COMMON_THROWEXCEPTION("There is no EZdown in an Seismic modelling")
    return (EZdown);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEXleft()
{
    COMMON_THROWEXCEPTION("There is no EXleft in an Seismic modelling")
    return (EXleft);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEYleft()
{
    COMMON_THROWEXCEPTION("There is no EYleft in an Seismic modelling")
    return (EYleft);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEZleft()
{
    COMMON_THROWEXCEPTION("There is no EZleft in an Seismic modelling")
    return (EZleft);
}

//! \brief Getter routine for EX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEXright()
{
    COMMON_THROWEXCEPTION("There is no EXright in an Seismic modelling")
    return (EXright);
}

//! \brief Getter routine for EY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEYright()
{
    COMMON_THROWEXCEPTION("There is no EYright in an Seismic modelling")
    return (EYright);
}

//! \brief Getter routine for EZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefEZright()
{
    COMMON_THROWEXCEPTION("There is no EZright in an Seismic modelling")
    return (EZright);
}

//! \brief Getter routine for RX Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRX()
{
    COMMON_THROWEXCEPTION("There is no RX in an Seismic modelling")
    return (RX);
}

//! \brief Getter routine for RY Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRY()
{
    COMMON_THROWEXCEPTION("There is no RY in an Seismic modelling")
    return (RY);
}

//! \brief Getter routine for RZ Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsSeismic<ValueType>::getRefRZ()
{
    COMMON_THROWEXCEPTION("There is no RZ in an Seismic modelling")
    return (RZ);
}

template class KITGPI::Wavefields::WavefieldsSeismic<float>;
template class KITGPI::Wavefields::WavefieldsSeismic<double>;
