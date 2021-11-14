#include "Wavefields2Dtmem.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dtmem<ValueType>::getContextPtr()
{
    return (HX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D tmem wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dtmem<ValueType>::FD2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "tmem";
    numDimension = 2;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    this->initWavefield(HX, ctx, dist);
    this->initWavefield(HY, ctx, dist);
    this->initWavefield(EZ, ctx, dist);
    this->initWavefield(EZup, ctx, dist);
    this->initWavefield(EZdown, ctx, dist);
    this->initWavefield(EZleft, ctx, dist);
    this->initWavefield(EZright, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD2Dtmem<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 3 Wavefields in 2D tmem modeling: HX, HY, EZ */
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
void KITGPI::Wavefields::FD2Dtmem<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 0:
        break;
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

        this->getCurl(derivatives, curl, model.getDielectricPermittivity());
        this->getDiv(derivatives, div, model.getVelocityEM());

        IO::writeVector(curl, fileName + ".curl." + timeStep, fileFormat);
        IO::writeVector(div, fileName + ".div." + timeStep, fileFormat);
        break;
    }
    case 4:
        IO::writeVector(EZ, fileName + ".EZ." + timeStep, fileFormat);
        IO::writeVector(EZup, fileName + ".EZ.up." + timeStep, fileFormat);
        IO::writeVector(EZdown, fileName + ".EZ.down." + timeStep, fileFormat);
        break;
    case 5:
        IO::writeVector(EZ, fileName + ".EZ." + timeStep, fileFormat);
        IO::writeVector(EZleft, fileName + ".EZ.left." + timeStep, fileFormat);
        IO::writeVector(EZright, fileName + ".EZ.right." + timeStep, fileFormat);
        break;
    default:
        COMMON_THROWEXCEPTION("Invalid snapType.")
    }
}

/*! \brief decompose wavefields to parts.
 \param decomposeType decomposeType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{    
    if (decomposeType > 0) {
        lama::DenseVector<ValueType> Poynting1;
        lama::DenseVector<ValueType> Poynting2;
        if (decomposeType == 1) {              
            Poynting1 = EZ;
            Poynting1 *= HX;
            
            auto const &Dyf = derivatives.getDyf();
            Poynting2 = Dyf * EZ;
            Poynting2 *= wavefieldsDerivative.getRefEZ();
            Poynting2 *= -1;
            
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            EZup.unaryOp(Poynting2, common::UnaryOp::SIGN);
            EZup -= 1;
            EZup.unaryOp(EZup, common::UnaryOp::ABS);
            EZup.unaryOp(EZup, common::UnaryOp::SIGN);
            EZup *= EZ;
            
            EZdown.unaryOp(Poynting2, common::UnaryOp::SIGN);
            EZdown += 1;
            EZdown.unaryOp(EZdown, common::UnaryOp::SIGN);
            EZdown *= EZ;
        } else if (decomposeType == 2) {            
            Poynting1 = EZ;
            Poynting1 *= HY;
            Poynting1 *= -1;
            
            auto const &Dxf = derivatives.getDxf();
            Poynting2 = Dxf * EZ;
            Poynting2 *= wavefieldsDerivative.getRefEZ();
            Poynting2 *= -1;
            
            Poynting2 *= Poynting1.maxNorm() / Poynting2.maxNorm();
            Poynting2 += Poynting1;
            
            EZleft.unaryOp(Poynting2, common::UnaryOp::SIGN);
            EZleft -= 1;
            EZleft.unaryOp(EZleft, common::UnaryOp::ABS);
            EZleft.unaryOp(EZleft, common::UnaryOp::SIGN);
            EZleft *= EZ;
                        
            EZright.unaryOp(Poynting2, common::UnaryOp::SIGN);
            EZright += 1;
            EZright.unaryOp(EZright, common::UnaryOp::SIGN);
            EZright *= EZ;  
        }
    }
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::resetWavefields()
{
    this->resetWavefield(HX);
    this->resetWavefield(HY);   
    this->resetWavefield(EZ);
}

/*! \brief Get numDimension (2)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD2Dtmem<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (tmem)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD2Dtmem<ValueType>::getEquationType() const
{
    return (equationType);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefHZ()
{
    COMMON_THROWEXCEPTION("There is no HZ wavefield in the 2D tmem case.")
    return (HZ);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefEX()
{
    COMMON_THROWEXCEPTION("There is no EX wavefield in the 2D tmem case.")
    return (EX);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefEY()
{
    COMMON_THROWEXCEPTION("There is no EY wavefield in the 2D tmem case.")
    return (EY);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefRX()
{
    COMMON_THROWEXCEPTION("There is no RX wavefield in the 2D tmem case.")
    return (RX);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefRY()
{
    COMMON_THROWEXCEPTION("There is no RY wavefield in the 2D tmem case.")
    return (RY);
}

//! \brief Not valid in the 2D tmem case
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> &KITGPI::Wavefields::FD2Dtmem<ValueType>::getRefRZ()
{
    COMMON_THROWEXCEPTION("There is no RZ wavefield in the 2D tmem case.")
    return (RZ);
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity)
{
    scai::lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    scai::lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();

    std::unique_ptr<lama::Vector<ValueType>> update_tmpPtr(HY.newVector());
    scai::lama::Vector<ValueType> &update_tmp = *update_tmpPtr;

    curl = Dyf * HX;
    update_tmp = Dxf * HY;
    curl -= update_tmp;

    update_tmp = scai::lama::sqrt(InverseDielectricPermittivity);
    curl *= update_tmp;
}

template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
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
KITGPI::Wavefields::FD2Dtmem<ValueType> KITGPI::Wavefields::FD2Dtmem<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD2Dtmem<ValueType> result;
    result.HX = this->HX * rhs;
    result.HY = this->HY * rhs;
    result.EZ = this->EZ * rhs;
    result.EZup = this->EZup * rhs;
    result.EZdown = this->EZdown * rhs;
    result.EZleft = this->EZleft * rhs;
    result.EZright = this->EZright * rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dtmem<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD2Dtmem<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dtmem<ValueType> KITGPI::Wavefields::FD2Dtmem<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dtmem<ValueType> KITGPI::Wavefields::FD2Dtmem<ValueType>::operator*(KITGPI::Wavefields::FD2Dtmem<ValueType> rhs)
{
    KITGPI::Wavefields::FD2Dtmem<ValueType> result;
    result.HX = this->HX * rhs.HX;
    result.HY = this->HY * rhs.HY;
    result.EZ = this->EZ * rhs.EZ;
    result.EZup = this->EZup * rhs.EZup;
    result.EZdown = this->EZdown * rhs.EZdown;
    result.EZleft = this->EZleft * rhs.EZleft;
    result.EZright = this->EZright * rhs.EZright;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD2Dtmem<ValueType> KITGPI::Wavefields::FD2Dtmem<ValueType>::operator*=(KITGPI::Wavefields::FD2Dtmem<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX = rhs.getRefHX();
    HY = rhs.getRefHY();
    EZ = rhs.getRefEZ();
    EZup = rhs.getRefEZup();
    EZdown = rhs.getRefEZdown();
    EZleft = rhs.getRefEZleft();
    EZright = rhs.getRefEZright();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX -= rhs.getRefHX();
    HY -= rhs.getRefHY();
    EZ -= rhs.getRefEZ();
    EZup -= rhs.getRefEZup();
    EZdown -= rhs.getRefEZdown();
    EZleft -= rhs.getRefEZleft();
    EZright -= rhs.getRefEZright();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX += rhs.getRefHX();
    HY += rhs.getRefHY();
    EZ += rhs.getRefEZ();
    EZup += rhs.getRefEZup();
    EZdown += rhs.getRefEZdown();
    EZleft += rhs.getRefEZleft();
    EZright += rhs.getRefEZright();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::timesAssign(ValueType rhs)
{
    HX *= rhs;
    HY *= rhs;
    EZ *= rhs;
    EZup *= rhs;
    EZdown *= rhs;
    EZleft *= rhs;
    EZright *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD2Dtmem<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    HX = lhs * rhs.getRefHX();
    HY = lhs * rhs.getRefHY();
    EZ = lhs * rhs.getRefEZ();
    EZup = lhs * rhs.getRefEZup();
    EZdown = lhs * rhs.getRefEZdown();
    EZleft = lhs * rhs.getRefEZleft();
    EZright = lhs * rhs.getRefEZright();
}

template class KITGPI::Wavefields::FD2Dtmem<float>;
template class KITGPI::Wavefields::FD2Dtmem<double>;
