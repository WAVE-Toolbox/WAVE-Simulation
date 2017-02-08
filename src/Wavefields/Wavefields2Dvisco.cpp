#include "Wavefields2Dvisco.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dvisco<ValueType>::getContextPtr()
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
KITGPI::Wavefields::FD2Dvisco<ValueType>::FD2Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
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
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D visco-elastic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D visco-elastic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D visco-elastic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D visco-elastic case.")
    return (Szz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D visco-elastic case.")
    return (Syz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D visco-elastic case.")
    return (Sxz);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getVZ()
{
    COMMON_THROWEXCEPTION("There is no VZ wavefield in the 2D visco-elastic case.")
    return (VZ);
}

//! \brief Not valid in the 2D visco-elastic case
template <typename ValueType>
lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dvisco<ValueType>::getP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D visco-elastic case.")
    return (P);
}

template class KITGPI::Wavefields::FD2Dvisco<double>;
template class KITGPI::Wavefields::FD2Dvisco<float>;
