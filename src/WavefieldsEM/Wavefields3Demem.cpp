#include "Wavefields3Demem.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD3Demem<ValueType>::getContextPtr()
{
    return (HX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D emem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Demem<ValueType>::FD3Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "emem";
    numDimension = 3;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(HX, ctx, dist);
    this->initWavefield(HY, ctx, dist);
    this->initWavefield(HZ, ctx, dist);
    this->initWavefield(EX, ctx, dist);
    this->initWavefield(EY, ctx, dist);
    this->initWavefield(EZ, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD3Demem<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 6 Wavefields in 3D emem modeling: EZ,EY,EX, HX, HY, HZ */
    IndexType numWavefields = 6;
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
void KITGPI::Wavefields::FD3Demem<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName;
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
        this->getDiv(derivatives, div, model.getVelocityEM());

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
void KITGPI::Wavefields::FD3Demem<ValueType>::decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{    
}
/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::resetWavefields()
{
    this->resetWavefield(HX);
    this->resetWavefield(HY);
    this->resetWavefield(HZ);
    this->resetWavefield(EX);
    this->resetWavefield(EY);
    this->resetWavefield(EZ);
}

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD3Demem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (emem)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD3Demem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 3D emem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Demem<ValueType>::getRefRX()
{
    COMMON_THROWEXCEPTION("There is no RX wavefield in the 3D emem case.")
    return (RX);
}

//! \brief Not valid in the 3D emem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Demem<ValueType>::getRefRY()
{
    COMMON_THROWEXCEPTION("There is no RY wavefield in the 3D emem case.")
    return (RY);
}

//! \brief Not valid in the 3D emem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD3Demem<ValueType>::getRefRZ()
{
    COMMON_THROWEXCEPTION("There is no RZ wavefield in the 3D emem case.")
    return (RZ);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity)
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
    curl *= InverseDielectricPermittivity;
    curl = scai::lama::sqrt(curl);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
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
KITGPI::Wavefields::FD3Demem<ValueType> KITGPI::Wavefields::FD3Demem<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD3Demem<ValueType> result;
    result.HX = this->HX * rhs;
    result.HY = this->HY * rhs;
    result.HZ = this->HZ * rhs;
    result.EX = this->EX * rhs;
    result.EY = this->EY * rhs;
    result.EZ = this->EZ * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Demem<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Demem<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Demem<ValueType> KITGPI::Wavefields::FD3Demem<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Demem<ValueType> KITGPI::Wavefields::FD3Demem<ValueType>::operator*(KITGPI::Wavefields::FD3Demem<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Demem<ValueType> result;
    result.HX = this->HX * rhs.HX;
    result.HY = this->HY * rhs.HY;
    result.HZ = this->HZ * rhs.HZ;
    result.EX = this->EX * rhs.EX;
    result.EY = this->EY * rhs.EY;
    result.EZ = this->EZ * rhs.EZ;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Demem<ValueType> KITGPI::Wavefields::FD3Demem<ValueType>::operator*=(KITGPI::Wavefields::FD3Demem<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX = rhs.getRefHX();
    HY = rhs.getRefHY();
    HZ = rhs.getRefHZ();
    EX = rhs.getRefEX();
    EY = rhs.getRefEY();
    EZ = rhs.getRefEZ();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX -= rhs.getRefHX();
    HY -= rhs.getRefHY();
    HZ -= rhs.getRefHZ();
    EX -= rhs.getRefEX();
    EY -= rhs.getRefEY();
    EZ -= rhs.getRefEZ();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX += rhs.getRefHX();
    HY += rhs.getRefHY();
    HZ += rhs.getRefHZ();
    EX += rhs.getRefEX();
    EY += rhs.getRefEY();
    EZ += rhs.getRefEZ();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::timesAssign(ValueType rhs)
{
    HX *= rhs;
    HY *= rhs;
    HZ *= rhs;
    EX *= rhs;
    EY *= rhs;
    EZ *= rhs;
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::timesAssign(scai::lama::DenseVector<ValueType> rhs)
{
    HX *= rhs;
    HY *= rhs;
    HZ *= rhs;
    EX *= rhs;
    EY *= rhs;
    EZ *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Demem<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX = lhs * rhs.getRefHX();
    HY = lhs * rhs.getRefHY();
    HZ = lhs * rhs.getRefHZ();
    EX = lhs * rhs.getRefEX();
    EY = lhs * rhs.getRefEY();
    EZ = lhs * rhs.getRefEZ();
}

template class KITGPI::Wavefields::FD3Demem<float>;
template class KITGPI::Wavefields::FD3Demem<double>;
