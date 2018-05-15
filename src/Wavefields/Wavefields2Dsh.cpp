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
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &/*derivatives*/, scai::lama::Vector<ValueType> const &/*SWaveModulus*/, scai::lama::Vector<ValueType> const &/*PWaveModulus*/, IndexType partitionedOut)
{
    std::string fileBaseName = baseName + type;
    
    switch(snapType){ 
      case 1:
	this->writeWavefield(VX, "VZ", fileBaseName, t, partitionedOut);
	break;
      case 2:
	this->writeWavefield(Syy, "Syz", fileBaseName, t, partitionedOut);
	this->writeWavefield(Sxy, "Sxz", fileBaseName, t, partitionedOut);
	break;
      case 3:
      {
          COMMON_THROWEXCEPTION("Not implemented in Wavefields2Dsh.");
      }
      default:
	COMMON_THROWEXCEPTION("Invalid snapType.");
    }
}


/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::resetWavefields()
{
    this->resetWavefield(VZ);
    this->resetWavefield(Sxz);
    this->resetWavefield(Syz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D sh case.")
    return (Sxx);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D sh case.")
    return (Syy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D sh case.")
    return (Szz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefSxy()
{
    COMMON_THROWEXCEPTION("There is no Sxy wavefield in the 2D sh case.")
    return (Sxy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefVX()
{
    COMMON_THROWEXCEPTION("There is no VX wavefield in the 2D sh case.")
    return (VX);
}
//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefVY()
{
    COMMON_THROWEXCEPTION("There is no VY wavefield in the 2D sh case.")
    return (VY);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D sh case.")
    return (P);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D sh case.")
    return (Rxx);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D sh case.")
    return (Ryy);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D sh case.")
    return (Rzz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D sh case.")
    return (Ryz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D sh case.")
    return (Rxz);
}

//! \brief Not valid in the 2D sh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dsh<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D sh case.")
    return (Rxy);
}


/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*(ValueType rhs)
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
KITGPI::Wavefields::FD2Dsh<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dsh<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dsh<ValueType> KITGPI::Wavefields::FD2Dsh<ValueType>::operator*=(ValueType rhs)
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

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ = rhs.getRefVZ();
    Syz = rhs.getRefSyz();
    Sxz = rhs.getRefSxz();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ -= rhs.getRefVZ();
    Syz -= rhs.getRefSyz();
    Sxz -= rhs.getRefSxz();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ += rhs.getRefVZ();
    Syz += rhs.getRefSyz();
    Sxz += rhs.getRefSxz();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::timesAssign(ValueType rhs)
{
    VZ *= rhs;
    Syz *= rhs;
    Sxz *= rhs;
}

template class KITGPI::Wavefields::FD2Dsh<float>;
template class KITGPI::Wavefields::FD2Dsh<double>;