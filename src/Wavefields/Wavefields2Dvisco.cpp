#include "Wavefields2Dvisco.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD2Dvisco<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoelastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType>::FD2Dvisco(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(Sxx, ctx, dist);
    this->initWavefield(Syy, ctx, dist);
    this->initWavefield(Sxy, ctx, dist);
    this->initWavefield(Rxx, ctx, dist);
    this->initWavefield(Ryy, ctx, dist);
    this->initWavefield(Rxy, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 /param type Type of the Seismogram
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VX, "VX", type, t);
    this->writeWavefield(VY, "VY", type, t);
    this->writeWavefield(Sxx, "Sxx", type, t);
    this->writeWavefield(Syy, "Syy", type, t);
    this->writeWavefield(Sxy, "Sxy", type, t);
    this->writeWavefield(Rxx, "Rxx", type, t);
    this->writeWavefield(Ryy, "Ryy", type, t);
    this->writeWavefield(Rxy, "Rxy", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Sxy);
    this->resetWavefield(Rxx);
    this->resetWavefield(Ryy);
    this->resetWavefield(Rxy);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D visco-elastic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D visco-elastic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D visco-elastic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D visco-elastic case.")
    return (Szz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D visco-elastic case.")
    return (Syz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D visco-elastic case.")
    return (Sxz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getVZ()
{
    COMMON_THROWEXCEPTION("There is no VZ wavefield in the 2D visco-elastic case.")
    return (VZ);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D visco-elastic case.")
    return (P);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType> KITGPI::Wavefields::FD2Dvisco<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Wavefields::FD2Dvisco<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.P = this->P * rhs;
    result.Sxx = this->Sxx * rhs;
    result.Syy = this->Syy * rhs;
    result.Sxy = this->Sxy * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Wavefields::FD2Dvisco<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType> KITGPI::Wavefields::FD2Dvisco<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType> KITGPI::Wavefields::FD2Dvisco<ValueType>::operator*(KITGPI::Wavefields::FD2Dvisco<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dvisco<ValueType> result;
    result.VX = this->VX * rhs.VX;
    result.VY = this->VY * rhs.VY;
    result.P = this->P * rhs.P;
    result.Sxx = this->Sxx * rhs.Sxx;
    result.Syy = this->Syy * rhs.Syy;
    result.Sxy = this->Sxy * rhs.Sxy;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType> KITGPI::Wavefields::FD2Dvisco<ValueType>::operator*=(KITGPI::Wavefields::FD2Dvisco<ValueType> rhs)
{
    return rhs * *this;
}

template class KITGPI::Wavefields::FD2Dvisco<double>;
template class KITGPI::Wavefields::FD2Dvisco<float>;
