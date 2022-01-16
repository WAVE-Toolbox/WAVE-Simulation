#include "Wavefields2Dacoustic.hpp"
#include "../IO/IO.hpp"

using namespace scai;

template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(P, ctx, dist);
    this->initWavefield(Pup, ctx, dist);
    this->initWavefield(Pdown, ctx, dist);
    this->initWavefield(Pleft, ctx, dist);
    this->initWavefield(Pright, ctx, dist);
    this->initWavefield(VXup, ctx, dist);
    this->initWavefield(VXdown, ctx, dist);
    this->initWavefield(VXleft, ctx, dist);
    this->initWavefield(VXright, ctx, dist);
    this->initWavefield(VYup, ctx, dist);
    this->initWavefield(VYdown, ctx, dist);
    this->initWavefield(VYleft, ctx, dist);
    this->initWavefield(VYright, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dacoustic<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 3 Wavefields in 2D acoustic modeling: P, Vx, Vy */
    IndexType numWavefields = 3;
    return (this->getMemoryUsage(dist, numWavefields));
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
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType>::FD2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "acoustic";
    numDimension = 2;
    init(ctx, dist, numRelaxationMechanisms_in);
}

/*! \brief override Methode tor write Wavefield Snapshot to file
 *
 *
 \param snapType Type of the wavefield snapshots 1=Velocities 2=pressure 3=div + curl
 \param baseName base name of the output file
 \param t Current Timestep
 \param derivatives derivatives object only used to output div/curl
 \param model model object only used to output div/curl
 \param fileFormat Output file format 
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const & /*derivatives*/, Modelparameter::Modelparameter<ValueType> const & /*model*/, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 0:
        break;
    case 1:
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(P, fileName + ".P." + timeStep, fileFormat);
        break;
    case 3:
        COMMON_THROWEXCEPTION("There is no curl or div of wavefield in the 2D acoustic case.")
        break;
    case 4:
        IO::writeVector(P, fileName + ".P." + timeStep, fileFormat);
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        IO::writeVector(Pup, fileName + ".P.up." + timeStep, fileFormat);
        IO::writeVector(Pdown, fileName + ".P.down." + timeStep, fileFormat);
        IO::writeVector(VXup, fileName + ".VX.up." + timeStep, fileFormat);
        IO::writeVector(VXdown, fileName + ".VX.down." + timeStep, fileFormat);
        IO::writeVector(VYup, fileName + ".VY.up." + timeStep, fileFormat);
        IO::writeVector(VYdown, fileName + ".VY.down." + timeStep, fileFormat);
        break;
    case 5:
        IO::writeVector(P, fileName + ".P." + timeStep, fileFormat);
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        IO::writeVector(Pleft, fileName + ".P.left." + timeStep, fileFormat);
        IO::writeVector(Pright, fileName + ".P.right." + timeStep, fileFormat);
        IO::writeVector(VXleft, fileName + ".VX.left." + timeStep, fileFormat);
        IO::writeVector(VXright, fileName + ".VX.right." + timeStep, fileFormat);
        IO::writeVector(VYleft, fileName + ".VY.left." + timeStep, fileFormat);
        IO::writeVector(VYright, fileName + ".VY.right." + timeStep, fileFormat);
        break;
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief decompose wavefields to parts by Poynting vector.
 \param decomposition decomposeWavefieldType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{    
    if (decomposition > 0) {
        lama::DenseVector<ValueType> Poynting1;
        lama::DenseVector<ValueType> Poynting2;
        if (decomposition == 1) {  
            // for all
            Poynting1 = -P;
            Poynting1 *= VY;
            
            // for P
            auto const &Dyf = derivatives.getDyf();
            Poynting2 = Dyf * P;
            Poynting2 *= wavefieldsDerivative.getRefP();
            Poynting2 *= -1;
                        
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            Pup.unaryOp(Poynting2, common::UnaryOp::SIGN);
            Pup -= 1;
            Pup.unaryOp(Pup, common::UnaryOp::ABS);
            Pup.unaryOp(Pup, common::UnaryOp::SIGN);
            Pup *= P;
            
            Pdown.unaryOp(Poynting2, common::UnaryOp::SIGN);
            Pdown += 1;
            Pdown.unaryOp(Pdown, common::UnaryOp::SIGN);
            Pdown *= P;
            
            // for VX
            Poynting2 = Dyf * VX;
            Poynting2 *= wavefieldsDerivative.getRefVX();
            Poynting2 *= -1;
                        
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            VXup.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VXup -= 1;
            VXup.unaryOp(VXup, common::UnaryOp::ABS);
            VXup.unaryOp(VXup, common::UnaryOp::SIGN);
            VXup *= VX;
            
            VXdown.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VXdown += 1;
            VXdown.unaryOp(VXdown, common::UnaryOp::SIGN);
            VXdown *= VX;
            
            // for VY
            Poynting2 = Dyf * VY;
            Poynting2 *= wavefieldsDerivative.getRefVY();
            Poynting2 *= -1;
                        
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            VYup.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VYup -= 1;
            VYup.unaryOp(VYup, common::UnaryOp::ABS);
            VYup.unaryOp(VYup, common::UnaryOp::SIGN);
            VYup *= VY;
            
            VYdown.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VYdown += 1;
            VYdown.unaryOp(VYdown, common::UnaryOp::SIGN);
            VYdown *= VY;
        } else if (decomposition == 2) {  
            // for all
            Poynting1 = -P;
            Poynting1 *= VX;

            // for P
            auto const &Dxf = derivatives.getDxf();
            Poynting2 = Dxf * P;
            Poynting2 *= wavefieldsDerivative.getRefP();
            Poynting2 *= -1;
            
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            Pleft.unaryOp(Poynting2, common::UnaryOp::SIGN);
            Pleft -= 1;
            Pleft.unaryOp(Pleft, common::UnaryOp::ABS);
            Pleft.unaryOp(Pleft, common::UnaryOp::SIGN);
            Pleft *= P;
                        
            Pright.unaryOp(Poynting2, common::UnaryOp::SIGN);
            Pright += 1;
            Pright.unaryOp(Pright, common::UnaryOp::SIGN);
            Pright *= P;  
            
            // for VX
            Poynting2 = Dxf * VX;
            Poynting2 *= wavefieldsDerivative.getRefVX();
            Poynting2 *= -1;
            
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            VXleft.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VXleft -= 1;
            VXleft.unaryOp(VXleft, common::UnaryOp::ABS);
            VXleft.unaryOp(VXleft, common::UnaryOp::SIGN);
            VXleft *= VX;
                        
            VXright.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VXright += 1;
            VXright.unaryOp(VXright, common::UnaryOp::SIGN);
            VXright *= VX; 
            
            // for VY
            Poynting2 = Dxf * VY;
            Poynting2 *= wavefieldsDerivative.getRefVY();
            Poynting2 *= -1;
            
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            VYleft.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VYleft -= 1;
            VYleft.unaryOp(VYleft, common::UnaryOp::ABS);
            VYleft.unaryOp(VYleft, common::UnaryOp::SIGN);
            VYleft *= VY;
                        
            VYright.unaryOp(Poynting2, common::UnaryOp::SIGN);
            VYright += 1;
            VYright.unaryOp(VYright, common::UnaryOp::SIGN);
            VYright *= VY;
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::resetWavefields()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(P);
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dacoustic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dacoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefVZ()
{
    COMMON_THROWEXCEPTION("There is no Vz wavefield in the 2D acoustic case.")
    return (VZ);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D acoustic case.")
    return (Sxx);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D acoustic case.")
    return (Syy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D acoustic case.")
    return (Szz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSyz()
{
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D acoustic case.")
    return (Syz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSxz()
{
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D acoustic case.")
    return (Sxz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefSxy()
{
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 2D acoustic case.")
    return (Sxy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D acoustic case.")
    return (Rxx);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D acoustic case.")
    return (Ryy);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D acoustic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D acoustic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D acoustic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D acoustic case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Dacoustic<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.P = this->P * rhs;
    result.Pup = this->Pup * rhs;
    result.Pdown = this->Pdown * rhs;
    result.Pleft = this->Pleft * rhs;
    result.Pright = this->Pright * rhs;
    result.VXup = this->VXup * rhs;
    result.VXdown = this->VXdown * rhs;
    result.VXleft = this->VXleft * rhs;
    result.VXright = this->VXright * rhs;
    result.VYup = this->VYup * rhs;
    result.VYdown = this->VYdown * rhs;
    result.VYleft = this->VYleft * rhs;
    result.VYright = this->VYright * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType> KITGPI::Wavefields::FD2Dacoustic<ValueType>::operator*=(ValueType rhs)
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
    result.Pup = this->Pup * rhs.Pup;
    result.Pdown = this->Pdown * rhs.Pdown;
    result.Pleft = this->Pleft * rhs.Pleft;
    result.Pright = this->Pright * rhs.Pright;
    result.VXup = this->VXup * rhs.VXup;
    result.VXdown = this->VXdown * rhs.VXdown;
    result.VXleft = this->VXleft * rhs.VXleft;
    result.VXright = this->VXright * rhs.VXright;
    result.VYup = this->VYup * rhs.VYup;
    result.VYdown = this->VYdown * rhs.VYdown;
    result.VYleft = this->VYleft * rhs.VYleft;
    result.VYright = this->VYright * rhs.VYright;
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

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = rhs.getRefVX();
    VY = rhs.getRefVY();
    P = rhs.getRefP();
    Pup = rhs.getRefPup();
    Pdown = rhs.getRefPdown();
    Pleft = rhs.getRefPleft();
    Pright = rhs.getRefPright();
    VXup = rhs.getRefVXup();
    VXdown = rhs.getRefVXdown();
    VXleft = rhs.getRefVXleft();
    VXright = rhs.getRefVXright();
    VYup = rhs.getRefVYup();
    VYdown = rhs.getRefVYdown();
    VYleft = rhs.getRefVYleft();
    VYright = rhs.getRefVYright();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX -= rhs.getRefVX();
    VY -= rhs.getRefVY();
    P -= rhs.getRefP();
    Pup -= rhs.getRefPup();
    Pdown -= rhs.getRefPdown();
    Pleft -= rhs.getRefPleft();
    Pright -= rhs.getRefPright();
    VXup -= rhs.getRefVXup();
    VXdown -= rhs.getRefVXdown();
    VXleft -= rhs.getRefVXleft();
    VXright -= rhs.getRefVXright();
    VYup -= rhs.getRefVYup();
    VYdown -= rhs.getRefVYdown();
    VYleft -= rhs.getRefVYleft();
    VYright -= rhs.getRefVYright();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX += rhs.getRefVX();
    VY += rhs.getRefVY();
    P += rhs.getRefP();
    Pup += rhs.getRefPup();
    Pdown += rhs.getRefPdown();
    Pleft += rhs.getRefPleft();
    Pright += rhs.getRefPright();
    VXup += rhs.getRefVXup();
    VXdown += rhs.getRefVXdown();
    VXleft += rhs.getRefVXleft();
    VXright += rhs.getRefVXright();
    VYup += rhs.getRefVYup();
    VYdown += rhs.getRefVYdown();
    VYleft += rhs.getRefVYleft();
    VYright += rhs.getRefVYright();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::timesAssign(ValueType rhs)
{
    VX *= rhs;
    VY *= rhs;
    P *= rhs;
    Pup *= rhs;
    Pdown *= rhs;
    Pleft *= rhs;
    Pright *= rhs;
    VXup *= rhs;
    VXdown *= rhs;
    VXleft *= rhs;
    VXright *= rhs;
    VYup *= rhs;
    VYdown *= rhs;
    VYleft *= rhs;
    VYright *= rhs;
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::timesAssign(scai::lama::DenseVector<ValueType> rhs)
{
    VX *= rhs;
    VY *= rhs;
    P *= rhs;
    Pup *= rhs;
    Pdown *= rhs;
    Pleft *= rhs;
    Pright *= rhs;
    VXup *= rhs;
    VXdown *= rhs;
    VXleft *= rhs;
    VXright *= rhs;
    VYup *= rhs;
    VYdown *= rhs;
    VYleft *= rhs;
    VYright *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = lhs * rhs.getRefVX();
    VY = lhs * rhs.getRefVY();
    P = lhs * rhs.getRefP();
    Pup = lhs * rhs.getRefPup();
    Pdown = lhs * rhs.getRefPdown();
    Pleft = lhs * rhs.getRefPleft();
    Pright = lhs * rhs.getRefPright();
    VXup = lhs * rhs.getRefVXup();
    VXdown = lhs * rhs.getRefVXdown();
    VXleft = lhs * rhs.getRefVXleft();
    VXright = lhs * rhs.getRefVXright();
    VYup = lhs * rhs.getRefVYup();
    VYdown = lhs * rhs.getRefVYdown();
    VYleft = lhs * rhs.getRefVYleft();
    VYright = lhs * rhs.getRefVYright();
}

template class KITGPI::Wavefields::FD2Dacoustic<double>;
template class KITGPI::Wavefields::FD2Dacoustic<float>;
