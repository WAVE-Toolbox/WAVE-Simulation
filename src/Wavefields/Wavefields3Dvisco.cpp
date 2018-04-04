#include "Wavefields3Dvisco.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD3Dvisco<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D viscoelastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType>::FD3Dvisco(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
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
    this->initWavefield(Rxx, ctx, dist);
    this->initWavefield(Ryy, ctx, dist);
    this->initWavefield(Rzz, ctx, dist);
    this->initWavefield(Ryz, ctx, dist);
    this->initWavefield(Rxz, ctx, dist);
    this->initWavefield(Rxy, ctx, dist);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param type Type of the Seismogram
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::write(IndexType snapType, std::string baseName,std::string type, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector const &SWaveModulus, scai::lama::Vector const &PWaveModulus, IndexType partitionedOut)
{
    std::string fileBaseName = baseName + type;
    
    switch(snapType){ 
      case 1:
	this->writeWavefield(VX, "VX", fileBaseName, t, partitionedOut);
	this->writeWavefield(VY, "VY", fileBaseName, t, partitionedOut);
	this->writeWavefield(VZ, "VZ", fileBaseName, t, partitionedOut);
	break;
      case 2:
	this->writeWavefield(Sxx, "Sxx", fileBaseName, t, partitionedOut);
	this->writeWavefield(Syy, "Syy", fileBaseName, t, partitionedOut);
	this->writeWavefield(Szz, "Szz", fileBaseName, t, partitionedOut);
	this->writeWavefield(Sxy, "Sxy", fileBaseName, t, partitionedOut);
	this->writeWavefield(Sxz, "Sxz", fileBaseName, t, partitionedOut);
	this->writeWavefield(Syz, "Syz", fileBaseName, t, partitionedOut);
	break;
      case 3:
      {
	common::unique_ptr<scai::lama::Vector> curl_Ptr(VX.newVector()); 
	scai::lama::Vector &curl = *curl_Ptr;
	common::unique_ptr<scai::lama::Vector> div_Ptr(VX.newVector()); 
	scai::lama::Vector &div = *div_Ptr;
	
	this->getCurl(derivatives,curl,SWaveModulus);
	this->getDiv(derivatives,div,PWaveModulus);
	
	scai::lama::DenseVector<ValueType>curlDense(curl);
	scai::lama::DenseVector<ValueType>divDense(div);
	this->writeWavefield(curlDense, "CURL", fileBaseName, t, partitionedOut);
	this->writeWavefield(divDense, "DIV", fileBaseName, t, partitionedOut);
      }
	break;
      default:
	COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief Wrapper Function to Write Snapshot of the Wavefield
 *
 *
 \param t Current Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::writeSnapshot(IndexType snapType, std::string baseName,IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector const &SWaveModulus, scai::lama::Vector const &PWaveModulus, IndexType partitionedOut)
{
    write(snapType, baseName, type, t, derivatives, SWaveModulus, PWaveModulus, partitionedOut);
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::resetWavefields()
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
    this->resetWavefield(Rxx);
    this->resetWavefield(Ryy);
    this->resetWavefield(Rzz);
    this->resetWavefield(Ryz);
    this->resetWavefield(Rxz);
    this->resetWavefield(Rxy);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector &curl, scai::lama::Vector const &SWaveModulus) 
{
    scai::lama::Matrix const &Dxb = derivatives.getDxb();
    scai::lama::Matrix const &Dyb = derivatives.getDzb();
    scai::lama::Matrix const &Dzb = derivatives.getDyb();
    
    common::unique_ptr<scai::lama::Vector> update_tmp1Ptr(VZ.newVector()); 
    scai::lama::Vector &update_tmp1 = *update_tmp1Ptr;  
    common::unique_ptr<scai::lama::Vector> update_tmp2Ptr(VZ.newVector()); 
    scai::lama::Vector &update_tmp2 = *update_tmp2Ptr;  
    
    common::unique_ptr<scai::lama::Vector> update_xPtr(VX.newVector()); 
    scai::lama::Vector &update_x = *update_xPtr;
    common::unique_ptr<scai::lama::Vector> update_yPtr(VY.newVector()); 
    scai::lama::Vector &update_y = *update_yPtr;   
    common::unique_ptr<scai::lama::Vector> update_zPtr(VZ.newVector()); 
    scai::lama::Vector &update_z = *update_zPtr;   
    
    //squared curl of velocity field
    update_tmp1 = Dyb * VZ;
    update_tmp2 = Dzb * VY;
    update_x = update_tmp1 - update_tmp2;
    update_x.powExp(2.0);
    curl = update_x;
    update_tmp1 = Dzb * VX;
    update_tmp2 = Dxb * VZ;
    update_y = update_tmp1 - update_tmp2;
    update_y.powExp(2.0);
    curl += update_y;
    update_tmp1 = Dxb * VY;
    update_tmp2 = Dyb * VX;
    update_z = update_tmp1 - update_tmp2;
    update_z.powExp(2.0);
    curl += update_z;
    
    // conversion to energy according to Dougherty and Stephen (PAGEOPH, 1988)
    curl *= SWaveModulus;
    curl.sqrt();
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector &div, lama::Vector const &PWaveModulus)
{
    scai::lama::Matrix const &Dxb = derivatives.getDxb();
    scai::lama::Matrix const &Dyb = derivatives.getDzb();
    scai::lama::Matrix const &Dzb = derivatives.getDyb();
    
    div = Dxb * VX;
    div += Dyb * VY;
    div += Dzb * VZ;
    
    div.powExp(2.0);
    div *= PWaveModulus;
    div.sqrt();
}

//! \brief Not valid in the 3D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dvisco<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D visco-elastic case.")
    return (P);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD3Dvisco<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.VZ = this->VZ * rhs;
    result.P = this->P * rhs;
    result.Sxx = this->Sxx * rhs;
    result.Syy = this->Syy * rhs;
    result.Szz = this->Szz * rhs;
    result.Sxy = this->Sxy * rhs;
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
KITGPI::Wavefields::FD3Dvisco<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Dvisco<ValueType> result;
    result.VX = this->VX * rhs.VX;
    result.VY = this->VY * rhs.VY;
    result.VZ = this->VZ * rhs.VZ;
    result.P = this->P * rhs.P;
    result.Sxx = this->Sxx * rhs.Sxx;
    result.Syy = this->Syy * rhs.Syy;
    result.Szz = this->Szz * rhs.Szz;
    result.Sxy = this->Sxy * rhs.Sxy;
    result.Sxz = this->Sxz * rhs.Sxz;
    result.Syz = this->Syz * rhs.Syz;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*=(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = rhs.getRefVX();
    VY = rhs.getRefVY();
    VZ = rhs.getRefVZ();
    Sxx = rhs.getRefSxx();
    Syy = rhs.getRefSyy();
    Szz = rhs.getRefSzz();
    Sxy = rhs.getRefSxy();
    Sxz = rhs.getRefSxz();
    Syz = rhs.getRefSyz();
    Rxx = rhs.getRefRxx();
    Ryy = rhs.getRefRyy();
    Rzz = rhs.getRefRzz();
    Rxy = rhs.getRefRxy();
    Rxz = rhs.getRefRxz();
    Ryz = rhs.getRefRyz();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX -= rhs.getRefVX();
    VY -= rhs.getRefVY();
    VZ -= rhs.getRefVZ();
    Sxx -= rhs.getRefSxx();
    Syy -= rhs.getRefSyy();
    Szz -= rhs.getRefSzz();
    Sxy -= rhs.getRefSxy();
    Sxz -= rhs.getRefSxz();
    Syz -= rhs.getRefSyz();
    Rxx -= rhs.getRefRxx();
    Ryy -= rhs.getRefRyy();
    Rzz -= rhs.getRefRzz();
    Rxy -= rhs.getRefRxy();
    Rxz -= rhs.getRefRxz();
    Ryz -= rhs.getRefRyz();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX += rhs.getRefVX();
    VY += rhs.getRefVY();
    VZ += rhs.getRefVZ();
    Sxx += rhs.getRefSxx();
    Syy += rhs.getRefSyy();
    Szz += rhs.getRefSzz();
    Sxy += rhs.getRefSxy();
    Sxz += rhs.getRefSxz();
    Syz += rhs.getRefSyz();
    Rxx += rhs.getRefRxx();
    Ryy += rhs.getRefRyy();
    Rzz += rhs.getRefRzz();
    Rxy += rhs.getRefRxy();
    Rxz += rhs.getRefRxz();
    Ryz += rhs.getRefRyz();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::timesAssign(ValueType rhs)
{
    VX *= rhs;
    VY *= rhs;
    VZ *= rhs;
    Sxx *= rhs;
    Syy *= rhs;
    Szz *= rhs;
    Sxy *= rhs;
    Sxz *= rhs;
    Syz *= rhs;
    Rxx *= rhs;
    Ryy *= rhs;
    Rzz *= rhs;
    Rxy *= rhs;
    Rxz *= rhs;
    Ryz *= rhs;
}

template class KITGPI::Wavefields::FD3Dvisco<float>;
template class KITGPI::Wavefields::FD3Dvisco<double>;
