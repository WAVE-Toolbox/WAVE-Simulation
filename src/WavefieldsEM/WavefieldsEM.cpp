#include "WavefieldsEM.hpp"
#include <math.h>

using namespace scai;

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
    auto read_EZ = EZ.getLocalValues();

    for (IndexType localIndex=0;localIndex<ownedIndeces.size();localIndex++) {
        if (read_EX.size()!=0){
            //SCAI_ASSERT_ERROR(isfinite(read_EX[localIndex]),"Infinite or NaN value in EX wavefield!" << localIndex);
            if (isfinite(read_EX[localIndex])==false){
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



/* Seismic */
//! \brief Getter routine for vX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVX()
{
    COMMON_THROWEXCEPTION("There is no VX in an EM modelling")
    return (VX);
}

//! \brief Getter routine for vY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVY()
{
    COMMON_THROWEXCEPTION("There is no VY in an EM modelling")
    return (VY);
}

//! \brief Getter routine for vZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVZ()
{
    COMMON_THROWEXCEPTION("There is no VZ in an EM modelling")
    return (VZ);
}

//! \brief Getter routine for Sxx wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx in an EM modelling")
    return (Sxx);
}

//! \brief Getter routine for Syy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy in an EM modelling")
    return (Syy);
}

//! \brief Getter routine for Szz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz in an EM modelling")
    return (Szz);
}

//! \brief Getter routine for Syz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz in an EM modelling")
    return (Syz);
}

//! \brief Getter routine for Sxz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz in an EM modelling")
    return (Sxz);
}

//! \brief Getter routine for Sxy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefSxy()
{
    COMMON_THROWEXCEPTION("There is no Sxy in an EM modelling")
    return (Sxy);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no P in an EM modelling")
    return (P);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefPup()
{
    COMMON_THROWEXCEPTION("There is no Pup in an EM modelling")
    return (Pup);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefPdown()
{
    COMMON_THROWEXCEPTION("There is no Pdown in an EM modelling")
    return (Pdown);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefPleft()
{
    COMMON_THROWEXCEPTION("There is no Pleft in an EM modelling")
    return (Pleft);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefPright()
{
    COMMON_THROWEXCEPTION("There is no Pright in an EM modelling")
    return (Pright);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVXup()
{
    COMMON_THROWEXCEPTION("There is no Vxup in an EM modelling")
    return (VXup);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVXdown()
{
    COMMON_THROWEXCEPTION("There is no VXdown in an EM modelling")
    return (VXdown);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVXleft()
{
    COMMON_THROWEXCEPTION("There is no VXleft in an EM modelling")
    return (VXleft);
}

//! \brief Getter routine for VX
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVXright()
{
    COMMON_THROWEXCEPTION("There is no VXright in an EM modelling")
    return (VXright);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVYup()
{
    COMMON_THROWEXCEPTION("There is no VYup in an EM modelling")
    return (VYup);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVYdown()
{
    COMMON_THROWEXCEPTION("There is no VYdown in an EM modelling")
    return (VYdown);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVYleft()
{
    COMMON_THROWEXCEPTION("There is no VYleft in an EM modelling")
    return (VYleft);
}

//! \brief Getter routine for VY
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVYright()
{
    COMMON_THROWEXCEPTION("There is no VYright in an EM modelling")
    return (VYright);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVZup()
{
    COMMON_THROWEXCEPTION("There is no VZup in an EM modelling")
    return (VZup);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVZdown()
{
    COMMON_THROWEXCEPTION("There is no VZdown in an EM modelling")
    return (VZdown);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVZleft()
{
    COMMON_THROWEXCEPTION("There is no VZleft in an EM modelling")
    return (VZleft);
}

//! \brief Getter routine for VZ
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefVZright()
{
    COMMON_THROWEXCEPTION("There is no VZright in an EM modelling")
    return (VZright);
}

//! \brief Getter routine for Rxx Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx in an EM modelling")
    return (Rxx);
}

//! \brief Getter routine for Ryy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy in an EM modelling")
    return (Ryy);
}

//! \brief Getter routine for Rzz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz in an EM modelling")
    return (Rzz);
}

//! \brief Getter routine for Ryz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz in an EM modelling")
    return (Ryz);
}

//! \brief Getter routine for Rxz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz in an EM modelling")
    return (Rxz);
}

//! \brief Getter routine for Rxy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::WavefieldsEM<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy in an EM modelling")
    return (Rxy);
}

template class KITGPI::Wavefields::WavefieldsEM<float>;
template class KITGPI::Wavefields::WavefieldsEM<double>;
