#include "Wavefields3Dacoustic.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD3Dacoustic<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D acoustic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType>::FD3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(P, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 /param type Type of the Seismogram
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::write(std::string type, IndexType t)
{
    this->writeWavefield(VX, "VX", type, t);
    this->writeWavefield(VY, "VY", type, t);
    this->writeWavefield(VZ, "VZ", type, t);
    this->writeWavefield(P, "P", type, t);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 /param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::writeSnapshot(IndexType t)
{
    write(type, t);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(P);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 3D acoustic case.")
    return (Sxx);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 3D acoustic case.")
    return (Syy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 3D acoustic case.")
    return (Szz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 3D acoustic case.")
    return (Syz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 3D acoustic case.")
    return (Sxz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxy()
{
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 3D acoustic case.")
    return (Sxy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 3D acoustic case.")
    return (Rxx);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 3D acoustic case.")
    return (Ryy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 3D acoustic case.")
    return (Rzz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 3D acoustic case.")
    return (Ryz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 3D acoustic case.")
    return (Rxz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 3D acoustic case.")
    return (Rxy);
}

template class KITGPI::Wavefields::FD3Dacoustic<float>;
template class KITGPI::Wavefields::FD3Dacoustic<double>;
