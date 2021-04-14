#include "Wavefields2Dviscotmem.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getContextPtr()
{
    return (HX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscotmem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType>::FD2Dviscotmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscotmem";
    numDimension = 2;
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(HX, ctx, dist);
    this->initWavefield(HY, ctx, dist);
    this->initWavefield(EZ, ctx, dist);
    this->initWavefield(RZ, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dviscotmem<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 4 Wavefields in 2D viscotmem modeling: HX, HY, EZ, RZ */
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
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName + type;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(HX, fileName + ".HX." + timeStep, fileFormat);
        IO::writeVector(HY, fileName + ".HY." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(EZ, fileName + ".EZ." + timeStep, fileFormat);
        break;
    case 3: {
        std::unique_ptr<lama::Vector<ValueType>> curl_Ptr(HX.newVector());
        scai::lama::Vector<ValueType> &curl = *curl_Ptr;
        std::unique_ptr<lama::Vector<ValueType>> div_Ptr(HX.newVector());
        scai::lama::Vector<ValueType> &div = *div_Ptr;

        this->getCurl(derivatives, curl, model.getDielectricPermittivityEM());
        this->getDiv(derivatives, div, model.getVelocityEM());

        IO::writeVector(curl, fileName + ".curl." + timeStep, fileFormat);
        IO::writeVector(div, fileName + ".div." + timeStep, fileFormat);
        break;
    }
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::resetWavefields()
{
    this->resetWavefield(HX);
    this->resetWavefield(HY);   
    this->resetWavefield(EZ); 
    this->resetWavefield(RZ);
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (tmem)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2D viscotmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getRefHZ()
{
    COMMON_THROWEXCEPTION("There is no HZ wavefield in the 2D viscotmem case.")
    return (HZ);
}

//! \brief Not valid in the 2D viscotmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getRefEX()
{
    COMMON_THROWEXCEPTION("There is no EX wavefield in the 2D viscotmem case.")
    return (EX);
}

//! \brief Not valid in the 2D viscotmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getRefEY()
{
    COMMON_THROWEXCEPTION("There is no EY wavefield in the 2D viscotmem case.")
    return (EY);
}

//! \brief Not valid in the 2D viscotmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getRefRX()
{
    COMMON_THROWEXCEPTION("There is no RX wavefield in the 2D viscotmem case.")
    return (RX);
}

//! \brief Not valid in the 2D viscotmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getRefRY()
{
    COMMON_THROWEXCEPTION("There is no RY wavefield in the 2D viscotmem case.")
    return (RY);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivityEM)
{
    scai::lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    scai::lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();

    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(HY.newVector());
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;

    curl = Dyf * HX;
    update_tmp = Dxf * HY;
    curl -= update_tmp;

    update_tmp = scai::lama::sqrt(InverseDielectricPermittivityEM);
    curl *= update_tmp;
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();

    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(HY.newVector());
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;

    div = Dxb * HX;
    div += Dyb * HY;

    update_tmp = scai::lama::sqrt(PWaveModulus);
    div *= update_tmp;
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType> KITGPI::Wavefields::FD2Dviscotmem<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Dviscotmem<ValueType> result;
    result.HX = this->HX * rhs;
    result.HY = this->HY * rhs;
    result.EZ = this->EZ * rhs;
    result.RZ = this->RZ * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dviscotmem<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType> KITGPI::Wavefields::FD2Dviscotmem<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType> KITGPI::Wavefields::FD2Dviscotmem<ValueType>::operator*(KITGPI::Wavefields::FD2Dviscotmem<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dviscotmem<ValueType> result;
    result.HX = this->HX * rhs.HX;
    result.HY = this->HY * rhs.HY;
    result.EZ = this->EZ * rhs.EZ;
    result.RZ = this->RZ * rhs.RZ;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dviscotmem<ValueType> KITGPI::Wavefields::FD2Dviscotmem<ValueType>::operator*=(KITGPI::Wavefields::FD2Dviscotmem<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX = rhs.getRefHX();
    HY = rhs.getRefHY();
    EZ = rhs.getRefEZ();
    RZ = rhs.getRefRZ();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX -= rhs.getRefHX();
    HY -= rhs.getRefHY();
    EZ -= rhs.getRefEZ();
    RZ -= rhs.getRefRZ();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX += rhs.getRefHX();
    HY += rhs.getRefHY();
    EZ += rhs.getRefEZ();
    RZ += rhs.getRefRZ();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::timesAssign(ValueType rhs)
{
    HX *= rhs;
    HY *= rhs;
    EZ *= rhs;
    RZ *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dviscotmem<ValueType>::applyWavefieldTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX = lhs * rhs.getRefHX();
    HY = lhs * rhs.getRefHY();
    EZ = lhs * rhs.getRefEZ();
    RZ = lhs * rhs.getRefRZ();
}

template class KITGPI::Wavefields::FD2Dviscotmem<float>;
template class KITGPI::Wavefields::FD2Dviscotmem<double>;
