#include "Wavefields3Delastic.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD3Delastic<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D elastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Delastic<ValueType>::FD3Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(Sxx, ctx, dist);
    this->initWavefield(Syy, ctx, dist);
    this->initWavefield(Szz, ctx, dist);
    this->initWavefield(Syz, ctx, dist);
    this->initWavefield(Sxz, ctx, dist);
    this->initWavefield(Sxy, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 /param type Type of the Seismogram
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Delastic<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VX, "VX", type, t);
    this->writeWavefield(VY, "VY", type, t);
    this->writeWavefield(VZ, "VZ", type, t);
    this->writeWavefield(Sxx, "Sxx", type, t);
    this->writeWavefield(Syy, "Syy", type, t);
    this->writeWavefield(Szz, "Szz", type, t);
    this->writeWavefield(Sxy, "Sxy", type, t);
    this->writeWavefield(Sxz, "Sxz", type, t);
    this->writeWavefield(Syz, "Syz", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Delastic<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Delastic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Szz);
    this->resetWavefield(Syz);
    this->resetWavefield(Sxz);
    this->resetWavefield(Sxy);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D elastic case.")
    return (P);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 3D elastic case.")
    return (Rxx);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 3D elastic case.")
    return (Ryy);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 3D elastic case.")
    return (Rzz);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 3D elastic case.")
    return (Ryz);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 3D elastic case.")
    return (Rxz);
}

//! \brief Not valid in the 3D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Delastic<ValueType>::getRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 3D elastic case.")
    return (Rxy);
}

template class KITGPI::Wavefields::FD3Delastic<float>;
template class KITGPI::Wavefields::FD3Delastic<double>;
