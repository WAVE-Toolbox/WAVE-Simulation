#include "Wavefields3Dacoustic.hpp"
#include "../IO/IO.hpp"

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
KITGPI::Wavefields::FD3Dacoustic<ValueType>::FD3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "acoustic";
    numDimension = 3;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(P, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD3Dacoustic<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 4 Wavefields in 3D acoustic modeling: P, Vx, Vy, Vz */
    IndexType numWavefields = 4;
    return (this->getMemoryUsage(dist, numWavefields));
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
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const & /*derivatives*/, Modelparameter::Modelparameter<ValueType> const & /*model*/, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        IO::writeVector(VZ, fileName + ".VZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(P, fileName + ".P." + timeStep, fileFormat);
        break;
    case 3:
        COMMON_THROWEXCEPTION("There is no curl or div of wavefield in the 3D acoustic case.")
        break;
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief decompose wavefields to parts.
 \param decomposition decomposeWavefieldType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{ 
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

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD3Dacoustic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD3Dacoustic<ValueType>::getEquationType() const
{
    return (equationType);
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
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 3D acoustic case.")
    return (Rxx);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 3D acoustic case.")
    return (Ryy);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 3D acoustic case.")
    return (Rzz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 3D acoustic case.")
    return (Ryz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 3D acoustic case.")
    return (Rxz);
}

//! \brief Not valid in the 3D acoustic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Dacoustic<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 3D acoustic case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*(ValueType rhs)
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
KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType> KITGPI::Wavefields::FD3Dacoustic<ValueType>::operator*=(ValueType rhs)
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

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::timesAssign(scai::lama::DenseVector<ValueType> rhs)
{
    VX *= rhs;
    VY *= rhs;
    VZ *= rhs;
    P *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = lhs * rhs.getRefVX();
    VY = lhs * rhs.getRefVY();
    VZ = lhs * rhs.getRefVZ();
    P = lhs * rhs.getRefP();
}

template class KITGPI::Wavefields::FD3Dacoustic<float>;
template class KITGPI::Wavefields::FD3Dacoustic<double>;
