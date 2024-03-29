#include "Wavefields2Delastic.hpp"
#include "../IO/IO.hpp"

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
KITGPI::Wavefields::FD2Delastic<ValueType>::FD2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "elastic";
    numDimension = 2;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(Sxx, ctx, dist);
    this->initWavefield(Syy, ctx, dist);
    this->initWavefield(Sxy, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Delastic<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 5 Wavefields in 2D elastic modeling: Sxx,Syy,Sxy, Vx, Vy */
    IndexType numWavefields = 5;
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
void KITGPI::Wavefields::FD2Delastic<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(Sxx, fileName + ".Sxx." + timeStep, fileFormat);
        IO::writeVector(Syy, fileName + ".Syy." + timeStep, fileFormat);
        IO::writeVector(Sxy, fileName + ".Sxy." + timeStep, fileFormat);
        break;
    case 3: {
        std::unique_ptr<lama::Vector<ValueType>> curl_Ptr(VX.newVector());
        scai::lama::Vector<ValueType> &curl = *curl_Ptr;
        std::unique_ptr<lama::Vector<ValueType>> div_Ptr(VX.newVector());
        scai::lama::Vector<ValueType> &div = *div_Ptr;

        this->getCurl(derivatives, curl, model.getSWaveModulus());
        this->getDiv(derivatives, div, model.getPWaveModulus());

        IO::writeVector(curl, fileName + ".curl." + timeStep, fileFormat);
        IO::writeVector(div, fileName + ".div." + timeStep, fileFormat);
        break;
    }
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
void KITGPI::Wavefields::FD2Delastic<ValueType>::decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{ 
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

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Delastic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (elastic)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Delastic<ValueType>::getEquationType() const
{
    return (equationType);
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
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxx()
{
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D elastic case.")
    return (Rxx);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRyy()
{
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D elastic case.")
    return (Ryy);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRzz()
{
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D elastic case.")
    return (Rzz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRyz()
{
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D elastic case.")
    return (Ryz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxz()
{
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D elastic case.")
    return (Rxz);
}

//! \brief Not valid in the 2D elastic case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Delastic<ValueType>::getRefRxy()
{
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D elastic case.")
    return (Rxy);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    scai::lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();

    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(VY.newVector());
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;

    curl = Dyf * VX;
    update_tmp = Dxf * VY;
    curl -= update_tmp;

    update_tmp = scai::lama::sqrt(SWaveModulus);
    curl *= update_tmp;
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();

    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(VY.newVector());
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;

    div = Dxb * VX;
    div += Dyb * VY;

    update_tmp = scai::lama::sqrt(PWaveModulus);
    div *= update_tmp;
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

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::timesAssign(scai::lama::DenseVector<ValueType> rhs)
{
    VX *= rhs;
    VY *= rhs;
    Sxx *= rhs;
    Syy *= rhs;
    Sxy *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = lhs * rhs.getRefVX();
    VY = lhs * rhs.getRefVY();
    Sxx = lhs * rhs.getRefSxx();
    Syy = lhs * rhs.getRefSyy();
    Sxy = lhs * rhs.getRefSxy();
}

template class KITGPI::Wavefields::FD2Delastic<float>;
template class KITGPI::Wavefields::FD2Delastic<double>;
