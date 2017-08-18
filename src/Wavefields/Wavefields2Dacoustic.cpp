#include "Wavefields2Dacoustic.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(P, ctx, dist);
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD2Dacoustic<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D acoustic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType>::FD2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 /param type Type of the Seismogram
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VX, "VX", type, t);
    this->writeWavefield(VY, "VY", type, t);
    this->writeWavefield(P, "P", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(P);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getVZ()
{
    COMMON_THROWEXCEPTION("There is no Vz wavefield in the 2D acoustic case.")
    return (VZ);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D acoustic case.")
    return (Sxx);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D acoustic case.")
    return (Syy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D acoustic case.")
    return (Szz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D acoustic case.")
    return (Syz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D acoustic case.")
    return (Sxz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxy()
{
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 2D acoustic case.")
    return (Sxy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D acoustic case.")
    return (Rxx);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D acoustic case.")
    return (Ryy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D acoustic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D acoustic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D acoustic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D acoustic case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Wavefields::FD2Dacoustic<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.P = this->P * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dacoustic<ValueType> result;
    result.VX = this->VX * rhs.VX;
    result.VY = this->VY * rhs.VY;
    result.P = this->P * rhs.P;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*=(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs)
{
    return rhs * *this;
}

template class KITGPI::Wavefields::FD2Dacoustic<double>;
template class KITGPI::Wavefields::FD2Dacoustic<float>;
