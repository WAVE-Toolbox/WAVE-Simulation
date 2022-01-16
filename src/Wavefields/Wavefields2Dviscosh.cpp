#include "Wavefields2Dviscosh.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dviscosh<ValueType>::getContextPtr()
{
    return (VZ.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscosh wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType>::FD2Dviscosh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "viscosh";
    numDimension = 2;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(Sxz, ctx, dist);
    this->initWavefield(Syz, ctx, dist);
    
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    scai::lama::DenseVector<ValueType> temp;
    this->initWavefield(temp, ctx, dist);
    Rxz.clear();
    Ryz.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz.push_back(temp);
        Ryz.push_back(temp);
    }
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dviscosh<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 3 Wavefields in 2D viscosh modeling: Sxz, Syz, Vz, Rxz, Ryz */
    IndexType numWavefields = 3 + 2 * numRelaxationMechanisms_in;
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
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const & /*derivatives*/, Modelparameter::Modelparameter<ValueType> const & /*model*/, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(VZ, fileName + ".VZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(Sxz, fileName + ".Sxz." + timeStep, fileFormat);
        IO::writeVector(Syz, fileName + ".Syz." + timeStep, fileFormat);
        for (int l=0; l<numRelaxationMechanisms; l++) {
            IO::writeVector(Rxz[l], fileName + ".Rxz" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Ryz[l], fileName + ".Ryz" + std::to_string(l+1) + "." + timeStep, fileFormat);
        }
        break;
    case 3: {
        COMMON_THROWEXCEPTION("Not implemented in Wavefields2Dviscosh.");
    }
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.");
    }
}

/*! \brief decompose wavefields to parts.
 \param decomposition decomposeWavefieldType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{ 
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::resetWavefields()
{
    this->resetWavefield(VZ);
    this->resetWavefield(Sxz);
    this->resetWavefield(Syz);
    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        this->resetWavefield(Rxz[l]);
        this->resetWavefield(Ryz[l]);
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dviscosh<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (sh)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dviscosh<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefSxx()
{
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D viscosh case.")
    return (Sxx);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefSyy()
{
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D viscosh case.")
    return (Syy);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefSzz()
{
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D viscosh case.")
    return (Szz);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefSxy()
{
    COMMON_THROWEXCEPTION("There is no Sxy wavefield in the 2D viscosh case.")
    return (Sxy);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefVX()
{
    COMMON_THROWEXCEPTION("There is no VX wavefield in the 2D viscosh case.")
    return (VX);
}
//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefVY()
{
    COMMON_THROWEXCEPTION("There is no VY wavefield in the 2D viscosh case.")
    return (VY);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D viscosh case.")
    return (P);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D viscosh case.")
    return (Rxx);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D viscosh case.")
    return (Ryy);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D viscosh case.")
    return (Rzz);
}

//! \brief Not valid in the 2D viscosh case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dviscosh<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D viscosh case.")
    return (Rxy);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType> KITGPI::Wavefields::FD2Dviscosh<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Dviscosh<ValueType> result;
    result.VZ = this->VZ * rhs;
    result.Sxz = this->Sxz * rhs;
    result.Syz = this->Syz * rhs;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        result.Rxz[l] = this->Rxz[l] * rhs;
        result.Ryz[l] = this->Ryz[l] * rhs;
    }
    
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dviscosh<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType> KITGPI::Wavefields::FD2Dviscosh<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType> KITGPI::Wavefields::FD2Dviscosh<ValueType>::operator*(KITGPI::Wavefields::FD2Dviscosh<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dviscosh<ValueType> result;
    result.VZ = this->VZ * rhs.VZ;
    result.Sxz = this->Sxz * rhs.Sxz;
    result.Syz = this->Syz * rhs.Syz;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        result.Rxz[l] = this->Rxz[l] * rhs.Rxz[l];
        result.Ryz[l] = this->Ryz[l] * rhs.Ryz[l];
    }
    
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscosh<ValueType> KITGPI::Wavefields::FD2Dviscosh<ValueType>::operator*=(KITGPI::Wavefields::FD2Dviscosh<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ = rhs.getRefVZ();
    Syz = rhs.getRefSyz();
    Sxz = rhs.getRefSxz();
    Rxz = rhs.getRefRxz();
    Ryz = rhs.getRefRyz();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ -= rhs.getRefVZ();
    Syz -= rhs.getRefSyz();
    Sxz -= rhs.getRefSxz();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz[l] -= rhs.getRefRxz()[l];
        Ryz[l] -= rhs.getRefRyz()[l];
    }
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ += rhs.getRefVZ();
    Syz += rhs.getRefSyz();
    Sxz += rhs.getRefSxz();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz[l] += rhs.getRefRxz()[l];
        Ryz[l] += rhs.getRefRyz()[l];
    }
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::timesAssign(ValueType rhs)
{
    VZ *= rhs;
    Syz *= rhs;
    Sxz *= rhs;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz[l] *= rhs;
        Ryz[l] *= rhs;
    }
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::timesAssign(scai::lama::DenseVector<ValueType> rhs)
{
    VZ *= rhs;
    Syz *= rhs;
    Sxz *= rhs;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz[l] *= rhs;
        Ryz[l] *= rhs;
    }
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscosh<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VZ = lhs * rhs.getRefVZ();
    Sxz = lhs * rhs.getRefSxz();
    Syz = lhs * rhs.getRefSyz();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxz[l] = lhs * rhs.getRefRxz()[l];
        Ryz[l] = lhs * rhs.getRefRyz()[l];
    }
}

template class KITGPI::Wavefields::FD2Dviscosh<float>;
template class KITGPI::Wavefields::FD2Dviscosh<double>;
