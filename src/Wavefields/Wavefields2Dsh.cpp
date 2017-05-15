#include "Wavefields2Dsh.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dsh<ValueType>::getContextPtr()
{
    return (VZ.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D sh wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType>::FD2Dsh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(Sxz, ctx, dist);
    this->initWavefield(Syz, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 /param type Type of the Seismogram
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VZ, "VZ", type, t);
    this->writeWavefield(Sxz, "Sxz", type, t);
    this->writeWavefield(Syz, "Syz", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::reset()
{
    this->resetWavefield(VZ);
    this->resetWavefield(Sxz);
    this->resetWavefield(Syz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D sh case.")
    return (Sxx);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D sh case.")
    return (Syy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D sh case.")
    return (Szz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getSxy()
{
    COMMON_THROWEXCEPTION("There is no Sxy wavefield in the 2D sh case.")
    return (Sxy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getVX()
{
    COMMON_THROWEXCEPTION("There is no VX wavefield in the 2D sh case.")
    return (VX);
}
//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getVY()
{
    COMMON_THROWEXCEPTION("There is no VY wavefield in the 2D sh case.")
    return (VY);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D sh case.")
    return (P);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D sh case.")
    return (Rxx);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D sh case.")
    return (Ryy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D sh case.")
    return (Rzz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D sh case.")
    return (Ryz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D sh case.")
    return (Rxz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D sh case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Wavefields::FD2Dsh<ValueType> result;
    result.VZ = this->VZ * rhs;
    result.Sxz = this->Sxz * rhs;
    result.Syz = this->Syz * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Wavefields::FD2Dsh<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*(KITGPI::Wavefields::FD2Dsh<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dsh<ValueType> result;
    result.VZ = this->VZ * rhs.VZ;
    result.Sxz = this->Sxz * rhs.Sxz;
    result.Syz = this->Syz * rhs.Syz;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*=(KITGPI::Wavefields::FD2Dsh<ValueType> rhs)
{
    return rhs * *this;
}

template class KITGPI::Wavefields::FD2Dsh<float>;
template class KITGPI::Wavefields::FD2Dsh<double>;
