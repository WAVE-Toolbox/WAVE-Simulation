#include "Wavefields2Dsh.hpp"
#include "../IO/IO.hpp"

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
    equationType = "sh";
    numDimension = 2;
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(Sxz, ctx, dist);
    this->initWavefield(Syz, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dsh<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 3 Wavefields in 2D sh modeling: Sxz, Syz, Vz */
    IndexType numWavefields = 3;
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
void KITGPI::Wavefields::FD2Dsh<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const & /*derivatives*/, Modelparameter::Modelparameter<ValueType> const & /*model*/, IndexType fileFormat)
{
    std::string fileName = baseName + type;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(VZ, fileName + ".VZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(Sxz, fileName + ".Sxz." + timeStep, fileFormat);
        IO::writeVector(Syz, fileName + ".Syz." + timeStep, fileFormat);
        break;
    case 3: {
        COMMON_THROWEXCEPTION("Not implemented in Wavefields2Dsh.");
    }
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.");
    }
}

/*! \brief decompose wavefields to parts.
 \param decomposeType decomposeType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{ 
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

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dsh<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (sh)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dsh<ValueType>::getEquationType() const
{
    return (equationType);
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

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dsh<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ = lhs * rhs.getRefVZ();
    Sxz = lhs * rhs.getRefSxz();
    Syz = lhs * rhs.getRefSyz();
}

template class KITGPI::Wavefields::FD2Dsh<float>;
template class KITGPI::Wavefields::FD2Dsh<double>;
