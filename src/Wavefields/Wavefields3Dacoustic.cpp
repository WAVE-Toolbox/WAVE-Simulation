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
 \param ctx Context
 \param dist Distribution
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
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::write(std::string baseName,std::string type, IndexType t, IndexType partitionedOut)
{
    std::string fileBaseName = baseName + type;
    this->writeWavefield(VX, "VX", fileBaseName, t, partitionedOut);
    this->writeWavefield(VY, "VY", fileBaseName, t, partitionedOut);
    this->writeWavefield(VZ, "VZ", fileBaseName, t, partitionedOut);
    this->writeWavefield(P, "P", fileBaseName, t, partitionedOut);
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::writeSnapshot(std::string baseName,IndexType t, IndexType partitionedOut)
{
    write(baseName, type, t, partitionedOut);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::resetWavefields()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(P);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 3D acoustic case.")
    return (Sxx);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 3D acoustic case.")
    return (Syy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 3D acoustic case.")
    return (Szz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 3D acoustic case.")
    return (Syz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 3D acoustic case.")
    return (Sxz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefSxy()
{
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 3D acoustic case.")
    return (Sxy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 3D acoustic case.")
    return (Rxx);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 3D acoustic case.")
    return (Ryy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 3D acoustic case.")
    return (Rzz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 3D acoustic case.")
    return (Ryz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 3D acoustic case.")
    return (Rxz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 3D acoustic case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Wavefields::FD3Dacoustic<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.VZ = this->VZ * rhs;
    result.P = this->P * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*(KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Dacoustic<ValueType> result;
    result.VX = this->VX * rhs.VX;
    result.VY = this->VY * rhs.VY;
    result.VZ = this->VZ * rhs.VZ;
    result.P = this->P * rhs.P;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*=(KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = rhs.getRefVX();
    VY = rhs.getRefVY();
    VZ = rhs.getRefVZ();
    P = rhs.getRefP();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX -= rhs.getRefVX();
    VY -= rhs.getRefVY();
    VZ -= rhs.getRefVZ();
    P -= rhs.getRefP();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX += rhs.getRefVX();
    VY += rhs.getRefVY();
    VZ += rhs.getRefVZ();
    P += rhs.getRefP();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::timesAssign(ValueType rhs)
{
    VX *= rhs;
    VY *= rhs;
    VZ *= rhs;
    P *= rhs;
}
template class KITGPI::Wavefields::FD3Dacoustic<float>;
template class KITGPI::Wavefields::FD3Dacoustic<double>;
