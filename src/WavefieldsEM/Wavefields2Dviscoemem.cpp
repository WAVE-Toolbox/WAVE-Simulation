#include "Wavefields2Dviscoemem.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getContextPtr()
{
    return (EX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoemem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType>::FD2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "viscoemem";
    numDimension = 2;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(HZ, ctx, dist);
    this->initWavefield(EX, ctx, dist);
    this->initWavefield(EY, ctx, dist);
    
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    scai::lama::DenseVector<ValueType> temp;
    this->initWavefield(temp, ctx, dist);
    RX.clear();
    RY.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        RX.push_back(temp);
        RY.push_back(temp);
    }
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dviscoemem<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 4 Wavefields in 2D viscoemem modeling: EX, RX, EY, RY, Hz */
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
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(HZ, fileName + ".HZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(EX, fileName + ".EX." + timeStep, fileFormat);
        IO::writeVector(EY, fileName + ".EY." + timeStep, fileFormat);
        for (int l=0; l<numRelaxationMechanisms; l++) {
            IO::writeVector(RX[l], fileName + ".RX" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(RY[l], fileName + ".RY" + std::to_string(l+1) + "." + timeStep, fileFormat);
        }
        break;
    case 3: {
        COMMON_THROWEXCEPTION("Not implemented in Wavefields2Dviscoemem.");
        break;
    }
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{    
}
/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::resetWavefields()
{
    this->resetWavefield(HZ);
    this->resetWavefield(EX);
    this->resetWavefield(EY);
    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        this->resetWavefield(RX[l]);
        this->resetWavefield(RY[l]);
    }
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2Dviscoemem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getRefHX()
{
    COMMON_THROWEXCEPTION("There is no HX wavefield in the 2Dviscoemem case.")
    return (HX);
}
//! \brief Not valid in the 2Dviscoemem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getRefHY()
{
    COMMON_THROWEXCEPTION("There is no HY wavefield in the 2Dviscoemem case.")
    return (HY);
}

//! \brief Not valid in the 2Dviscoemem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getRefEZ()
{
    COMMON_THROWEXCEPTION("There is no EZ wavefield in the 2Dviscoemem case.")
    return (EZ);
}

//! \brief Not valid in the 2Dviscoemem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dviscoemem<ValueType>::getRefRZ()
{
    COMMON_THROWEXCEPTION("There is no RZ wavefield in the 2Dviscoemem case.")
    return (RZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType> KITGPI::Wavefields::FD2Dviscoemem<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Dviscoemem<ValueType> result;
    result.HZ = this->HZ * rhs;
    result.EX = this->EX * rhs;
    result.EY = this->EY * rhs;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        result.RX[l] = this->RX[l] * rhs;
        result.RY[l] = this->RY[l] * rhs;
    }
    
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType> KITGPI::Wavefields::FD2Dviscoemem<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType> KITGPI::Wavefields::FD2Dviscoemem<ValueType>::operator*(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dviscoemem<ValueType> result;
    result.HZ = this->HZ * rhs.HZ;
    result.EX = this->EX * rhs.EX;
    result.EY = this->EY * rhs.EY;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        result.RX[l] = this->RX[l] * rhs.RX[l];
        result.RY[l] = this->RY[l] * rhs.RY[l];
    }
    
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscoemem<ValueType> KITGPI::Wavefields::FD2Dviscoemem<ValueType>::operator*=(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HZ = rhs.getRefHZ();
    EX = rhs.getRefEX();
    EY = rhs.getRefEY();
    RX = rhs.getRefRX();
    RY = rhs.getRefRY();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HZ -= rhs.getRefHZ();
    EX -= rhs.getRefEX();
    EY -= rhs.getRefEY();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        RX[l] -= rhs.getRefRX()[l];
        RY[l] -= rhs.getRefRY()[l];
    }
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HZ += rhs.getRefHZ();
    EX += rhs.getRefEX();
    EY += rhs.getRefEY();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        RX[l] += rhs.getRefRX()[l];
        RY[l] += rhs.getRefRY()[l];
    }
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::timesAssign(ValueType rhs)
{
    HZ *= rhs;
    EX *= rhs;
    EY *= rhs;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        RX[l] *= rhs;
        RY[l] *= rhs;
    }
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscoemem<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HZ = lhs * rhs.getRefHZ();
    EX = lhs * rhs.getRefEX();
    EY = lhs * rhs.getRefEY();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        RX[l] = lhs * rhs.getRefRX()[l];
        RY[l] = lhs * rhs.getRefRY()[l];
    }
}

template class KITGPI::Wavefields::FD2Dviscoemem<double>;
template class KITGPI::Wavefields::FD2Dviscoemem<float>;
