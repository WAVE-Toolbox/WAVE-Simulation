#include "Wavefields3Dviscoemem.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD3Dviscoemem<ValueType>::getContextPtr()
{
    return (HX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D viscoemem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType>::FD3Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscoemem";
    numDimension = 3;
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(HX, ctx, dist);
    this->initWavefield(HY, ctx, dist);
    this->initWavefield(HZ, ctx, dist);
    this->initWavefield(EX, ctx, dist);
    this->initWavefield(EY, ctx, dist);
    this->initWavefield(EZ, ctx, dist);
    this->initWavefield(RX, ctx, dist);
    this->initWavefield(RY, ctx, dist);
    this->initWavefield(RZ, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD3Dviscoemem<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 9 Wavefields in 3D acoustic modeling: EZ,EY,EX,RZ,RY,RX, Hx, Hy, Hz */
    IndexType numWavefields = 9;
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
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName + type;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(HX, fileName + ".HX." + timeStep, fileFormat);
        IO::writeVector(HY, fileName + ".HY." + timeStep, fileFormat);
        IO::writeVector(HZ, fileName + ".HZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(EX, fileName + ".EX." + timeStep, fileFormat);
        IO::writeVector(EY, fileName + ".EY." + timeStep, fileFormat);
        IO::writeVector(EZ, fileName + ".EZ." + timeStep, fileFormat);
        break;
    case 3: {
        std::unique_ptr<lama::Vector<ValueType>> curl_Ptr(HX.newVector());
        scai::lama::Vector<ValueType> &curl = *curl_Ptr;
        std::unique_ptr<lama::Vector<ValueType>> div_Ptr(HX.newVector());
        scai::lama::Vector<ValueType> &div = *div_Ptr;

        this->getCurl(derivatives, curl, model.getDielectricPermittivity());
        this->getDiv(derivatives, div, model.getElectricConductivity());

        IO::writeVector(curl, fileName + ".CURL." + timeStep, fileFormat);
        IO::writeVector(div, fileName + ".DIV." + timeStep, fileFormat);
    } break;
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::decompose(IndexType decomposeType, KITGPI::Wavefields::WavefieldsEM<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{    
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::resetWavefields()
{
    this->resetWavefield(HX);
    this->resetWavefield(HY);
    this->resetWavefield(HZ);
    this->resetWavefield(EX);
    this->resetWavefield(EY);
    this->resetWavefield(EZ);
    this->resetWavefield(RX);
    this->resetWavefield(RY);
    this->resetWavefield(RZ);
}

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD3Dviscoemem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD3Dviscoemem<ValueType>::getEquationType() const
{
    return (equationType);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    scai::lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
    scai::lama::Matrix<ValueType> const &Dzf = derivatives.getDzf();

    std::unique_ptr<lama::Vector<ValueType>> update_Ptr(HX.newVector());
    scai::lama::Vector<ValueType> &update = *update_Ptr;

    //squared curl of velocity field
    update = Dyf * HZ;
    update -= Dzf * HY;
    update = scai::lama::pow(update, 2.0);
    curl = update;
    update = Dzf * HX;
    update -= Dxf * HZ;
    update = scai::lama::pow(update, 2.0);
    curl += update;
    update = Dxf * HY;
    update -= Dyf * HX;
    update = scai::lama::pow(update, 2.0);
    curl += update;

    // conversion to energy according to Dougherty and Stephen (PAGEOPH, 1988)
    curl *= SWaveModulus;
    curl = scai::lama::sqrt(curl);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    scai::lama::Matrix<ValueType> const &Dzb = derivatives.getDzb();

    div = Dxb * HX;
    div += Dyb * HY;
    div += Dzb * HZ;

    div = scai::lama::pow(div, 2.0);
    div *= PWaveModulus;
    div = scai::lama::sqrt(div);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType> KITGPI::Wavefields::FD3Dviscoemem<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD3Dviscoemem<ValueType> result;
    result.HX = this->HX * rhs;
    result.HY = this->HY * rhs;
    result.HZ = this->HZ * rhs;
    result.EX = this->EX * rhs;
    result.EY = this->EY * rhs;
    result.EZ = this->EZ * rhs;
    result.RX = this->RX * rhs;
    result.RY = this->RY * rhs;
    result.RZ = this->RZ * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Dviscoemem<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType> KITGPI::Wavefields::FD3Dviscoemem<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType> KITGPI::Wavefields::FD3Dviscoemem<ValueType>::operator*(KITGPI::Wavefields::FD3Dviscoemem<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Dviscoemem<ValueType> result;
    result.HX = this->HX * rhs.HX;
    result.HY = this->HY * rhs.HY;
    result.HZ = this->HZ * rhs.HZ;
    result.EX = this->EX * rhs.EX;
    result.EY = this->EY * rhs.EY;
    result.EZ = this->EZ * rhs.EZ;
    result.RX = this->RX * rhs.RX;
    result.RY = this->RY * rhs.RY;
    result.RZ = this->RZ * rhs.RZ;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoemem<ValueType> KITGPI::Wavefields::FD3Dviscoemem<ValueType>::operator*=(KITGPI::Wavefields::FD3Dviscoemem<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX = rhs.getRefHX();
    HY = rhs.getRefHY();
    HZ = rhs.getRefHZ();
    EX = rhs.getRefEX();
    EY = rhs.getRefEY();
    EZ = rhs.getRefEZ();
    RX = rhs.getRefRX();
    RY = rhs.getRefRY();
    RZ = rhs.getRefRZ();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX -= rhs.getRefHX();
    HY -= rhs.getRefHY();
    HZ -= rhs.getRefHZ();
    EX -= rhs.getRefEX();
    EY -= rhs.getRefEY();
    EZ -= rhs.getRefEZ();
    RX -= rhs.getRefRX();
    RY -= rhs.getRefRY();
    RZ -= rhs.getRefRZ();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX += rhs.getRefHX();
    HY += rhs.getRefHY();
    HZ += rhs.getRefHZ();
    EX += rhs.getRefEX();
    EY += rhs.getRefEY();
    EZ += rhs.getRefEZ();
    RX += rhs.getRefRX();
    RY += rhs.getRefRY();
    RZ += rhs.getRefRZ();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::timesAssign(ValueType rhs)
{
    HX *= rhs;
    HY *= rhs;
    HZ *= rhs;
    EX *= rhs;
    EY *= rhs;
    EZ *= rhs;
    RX *= rhs;
    RY *= rhs;
    RZ *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoemem<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs)
{
    HX = lhs * rhs.getRefHX();
    HY = lhs * rhs.getRefHY();
    HZ = lhs * rhs.getRefHZ();
    EX = lhs * rhs.getRefEX();
    EY = lhs * rhs.getRefEY();
    EZ = lhs * rhs.getRefEZ();
    RX = lhs * rhs.getRefRX();
    RY = lhs * rhs.getRefRY();
    RZ = lhs * rhs.getRefRZ();
}

template class KITGPI::Wavefields::FD3Dviscoemem<float>;
template class KITGPI::Wavefields::FD3Dviscoemem<double>;
