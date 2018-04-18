#include "Wavefields2Delastic.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Delastic<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D elastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Delastic<ValueType>::FD2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(Sxx, ctx, dist);
    this->initWavefield(Syy, ctx, dist);
    this->initWavefield(Sxy, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> const &SWaveModulus, scai::lama::Vector<ValueType> const &PWaveModulus, IndexType partitionedOut)
{
    std::string fileBaseName = baseName + type;
    
    switch(snapType){ 
      case 1:
	this->writeWavefield(VX, "VX", fileBaseName, t, partitionedOut);
	this->writeWavefield(VY, "VY", fileBaseName, t, partitionedOut);
	break;
      case 2:
	this->writeWavefield(Sxx, "Sxx", fileBaseName, t, partitionedOut);
	this->writeWavefield(Syy, "Syy", fileBaseName, t, partitionedOut);
	this->writeWavefield(Sxy, "Sxy", fileBaseName, t, partitionedOut);
	break;
      case 3:
      {
	std::unique_ptr<lama::Vector<ValueType>> curl_Ptr(VX.newVector()); 
	scai::lama::Vector<ValueType> &curl = *curl_Ptr;
	std::unique_ptr<lama::Vector<ValueType>> div_Ptr(VX.newVector()); 
	scai::lama::Vector<ValueType> &div = *div_Ptr;
	
	this->getCurl(derivatives,curl,SWaveModulus);
	this->getDiv(derivatives,div,PWaveModulus);
	
	this->writeWavefield(curl, "CURL", fileBaseName, t, partitionedOut);
	this->writeWavefield(div, "DIV", fileBaseName, t, partitionedOut);
	break;
      }
      default:
	COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::resetWavefields()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Sxy);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D elastic case.")
    return (Sxz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D elastic case.")
    return (Syz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D elastic case.")
    return (Szz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefVZ()
{
    COMMON_THROWEXCEPTION("There is no VZ wavefield in the 2D elastic case.")
    return (VZ);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D elastic case.")
    return (P);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D elastic case.")
    return (Rxx);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D elastic case.")
    return (Ryy);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D elastic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D elastic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D elastic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D elastic case.")
    return (Rxy);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    
    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(VX.newVector()); 
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;  
    
    curl = Dxb * VY;
    update_tmp = Dyb * VX;
    curl -= update_tmp;
    
    curl = scai::lama::pow(curl,2.0);
    curl *= SWaveModulus;
    curl = scai::lama::sqrt(curl);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType>const &PWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    
    div = Dxb * VX;
    div += Dyb * VY;
    
    div = scai::lama::pow(div,2.0);
    div *= PWaveModulus;
    div = scai::lama::sqrt(div);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Delastic<ValueType> KITGPI::Wavefields::FD2Delastic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Delastic<ValueType> result;
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
KITGPI::Wavefields::FD2Delastic<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Delastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Delastic<ValueType> KITGPI::Wavefields::FD2Delastic<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Delastic<ValueType> KITGPI::Wavefields::FD2Delastic<ValueType>::operator*(KITGPI::Wavefields::FD2Delastic<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Delastic<ValueType> result;
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
KITGPI::Wavefields::FD2Delastic<ValueType> KITGPI::Wavefields::FD2Delastic<ValueType>::operator*=(KITGPI::Wavefields::FD2Delastic<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = rhs.getRefVX();
    VY = rhs.getRefVY();
    Sxx = rhs.getRefSxx();
    Syy = rhs.getRefSyy();
    Sxy = rhs.getRefSxy();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX -= rhs.getRefVX();
    VY -= rhs.getRefVY();
    Sxx -= rhs.getRefSxx();
    Syy -= rhs.getRefSyy();
    Sxy -= rhs.getRefSxy();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX += rhs.getRefVX();
    VY += rhs.getRefVY();
    Sxx += rhs.getRefSxx();
    Syy += rhs.getRefSyy();
    Sxy += rhs.getRefSxy();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::timesAssign(ValueType rhs)
{
    VX *= rhs;
    VY *= rhs;
    Sxx *= rhs;
    Syy *= rhs;
    Sxy *= rhs;
}
template class KITGPI::Wavefields::FD2Delastic<float>;
template class KITGPI::Wavefields::FD2Delastic<double>;
