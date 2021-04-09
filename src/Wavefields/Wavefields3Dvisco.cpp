#include "Wavefields3Dvisco.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD3Dvisco<ValueType>::getContextPtr()
{
    return (VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D viscoelastic wavefields
 *
 \param ctx Context
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType>::FD3Dvisco(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscoelastic";
    numDimension = 3;
    init(ctx, dist);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    this->initWavefield(VX, ctx, dist);
    this->initWavefield(VY, ctx, dist);
    this->initWavefield(VZ, ctx, dist);
    this->initWavefield(Sxx, ctx, dist);
    this->initWavefield(Syy, ctx, dist);
    this->initWavefield(Szz, ctx, dist);
    this->initWavefield(Syz, ctx, dist);
    this->initWavefield(Sxz, ctx, dist);
    this->initWavefield(Sxy, ctx, dist);
    this->initWavefield(Rxx, ctx, dist);
    this->initWavefield(Ryy, ctx, dist);
    this->initWavefield(Rzz, ctx, dist);
    this->initWavefield(Ryz, ctx, dist);
    this->initWavefield(Rxz, ctx, dist);
    this->initWavefield(Rxy, ctx, dist);
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD3Dvisco<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 15 Wavefields in 2D acoustic modeling: Sxx,Syy,Szz,Sxy,Sxz,Syz,Rxx,Ryy,Rzz,Rxy,Rxz,Ryz, Vx, Vy, Vz */
    IndexType numWavefields = 15;
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
void KITGPI::Wavefields::FD3Dvisco<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
{
    std::string fileName = baseName + type;
    std::string timeStep = std::to_string(static_cast<long long>(t));

    switch (snapType) {
    case 1:
        IO::writeVector(VX, fileName + ".VX." + timeStep, fileFormat);
        IO::writeVector(VY, fileName + ".VY." + timeStep, fileFormat);
        IO::writeVector(VZ, fileName + ".VZ." + timeStep, fileFormat);
        break;
    case 2:
        IO::writeVector(Sxx, fileName + ".Sxx." + timeStep, fileFormat);
        IO::writeVector(Syy, fileName + ".Syy." + timeStep, fileFormat);
        IO::writeVector(Szz, fileName + ".Szz." + timeStep, fileFormat);
        IO::writeVector(Sxy, fileName + ".Sxy." + timeStep, fileFormat);
        IO::writeVector(Sxz, fileName + ".Sxz." + timeStep, fileFormat);
        IO::writeVector(Syz, fileName + ".Syz." + timeStep, fileFormat);
        break;
    case 3: {
        std::unique_ptr<lama::Vector<ValueType>> curl_Ptr(VX.newVector());
        scai::lama::Vector<ValueType> &curl = *curl_Ptr;
        std::unique_ptr<lama::Vector<ValueType>> div_Ptr(VX.newVector());
        scai::lama::Vector<ValueType> &div = *div_Ptr;

        this->getCurl(derivatives, curl, model.getSWaveModulus());
        this->getDiv(derivatives, div, model.getPWaveModulus());

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
void KITGPI::Wavefields::FD3Dvisco<ValueType>::resetWavefields()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Szz);
    this->resetWavefield(Syz);
    this->resetWavefield(Sxz);
    this->resetWavefield(Sxy);
    this->resetWavefield(Rxx);
    this->resetWavefield(Ryy);
    this->resetWavefield(Rzz);
    this->resetWavefield(Ryz);
    this->resetWavefield(Rxz);
    this->resetWavefield(Rxy);
}

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD3Dvisco<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscoelastic)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD3Dvisco<ValueType>::getEquationType() const
{
    return (equationType);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    scai::lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
    scai::lama::Matrix<ValueType> const &Dzf = derivatives.getDzf();

    std::unique_ptr<lama::Vector<ValueType>> update_Ptr(VX.newVector());
    scai::lama::Vector<ValueType> &update = *update_Ptr;

    //squared curl of velocity field
    update = Dyf * VZ;
    update -= Dzf * VY;
    update = scai::lama::pow(update, 2.0);
    curl = update;
    update = Dzf * VX;
    update -= Dxf * VZ;
    update = scai::lama::pow(update, 2.0);
    curl += update;
    update = Dxf * VY;
    update -= Dyf * VX;
    update = scai::lama::pow(update, 2.0);
    curl += update;

    // conversion to energy according to Dougherty and Stephen (PAGEOPH, 1988)
    curl *= SWaveModulus;
    curl = scai::lama::sqrt(curl);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
{
    scai::lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    scai::lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    scai::lama::Matrix<ValueType> const &Dzb = derivatives.getDzb();

    div = Dxb * VX;
    div += Dyb * VY;
    div += Dzb * VZ;

    div = scai::lama::pow(div, 2.0);
    div *= PWaveModulus;
    div = scai::lama::sqrt(div);
}

//! \brief Not valid in the 3D visco-elastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dvisco<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D visco-elastic case.")
    return (P);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD3Dvisco<ValueType> result;
    result.VX = this->VX * rhs;
    result.VY = this->VY * rhs;
    result.VZ = this->VZ * rhs;
    result.P = this->P * rhs;
    result.Sxx = this->Sxx * rhs;
    result.Syy = this->Syy * rhs;
    result.Szz = this->Szz * rhs;
    result.Sxy = this->Sxy * rhs;
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
KITGPI::Wavefields::FD3Dvisco<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Dvisco<ValueType> result;
    result.VX = this->VX * rhs.VX;
    result.VY = this->VY * rhs.VY;
    result.VZ = this->VZ * rhs.VZ;
    result.P = this->P * rhs.P;
    result.Sxx = this->Sxx * rhs.Sxx;
    result.Syy = this->Syy * rhs.Syy;
    result.Szz = this->Szz * rhs.Szz;
    result.Sxy = this->Sxy * rhs.Sxy;
    result.Sxz = this->Sxz * rhs.Sxz;
    result.Syz = this->Syz * rhs.Syz;
    return result;
}

/*! \brief Overloading *= Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType> KITGPI::Wavefields::FD3Dvisco<ValueType>::operator*=(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = rhs.getRefVX();
    VY = rhs.getRefVY();
    VZ = rhs.getRefVZ();
    Sxx = rhs.getRefSxx();
    Syy = rhs.getRefSyy();
    Szz = rhs.getRefSzz();
    Sxy = rhs.getRefSxy();
    Sxz = rhs.getRefSxz();
    Syz = rhs.getRefSyz();
    Rxx = rhs.getRefRxx();
    Ryy = rhs.getRefRyy();
    Rzz = rhs.getRefRzz();
    Rxy = rhs.getRefRxy();
    Rxz = rhs.getRefRxz();
    Ryz = rhs.getRefRyz();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is substracted.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX -= rhs.getRefVX();
    VY -= rhs.getRefVY();
    VZ -= rhs.getRefVZ();
    Sxx -= rhs.getRefSxx();
    Syy -= rhs.getRefSyy();
    Szz -= rhs.getRefSzz();
    Sxy -= rhs.getRefSxy();
    Sxz -= rhs.getRefSxz();
    Syz -= rhs.getRefSyz();
    Rxx -= rhs.getRefRxx();
    Ryy -= rhs.getRefRyy();
    Rzz -= rhs.getRefRzz();
    Rxy -= rhs.getRefRxy();
    Rxz -= rhs.getRefRxz();
    Ryz -= rhs.getRefRyz();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX += rhs.getRefVX();
    VY += rhs.getRefVY();
    VZ += rhs.getRefVZ();
    Sxx += rhs.getRefSxx();
    Syy += rhs.getRefSyy();
    Szz += rhs.getRefSzz();
    Sxy += rhs.getRefSxy();
    Sxz += rhs.getRefSxz();
    Syz += rhs.getRefSyz();
    Rxx += rhs.getRefRxx();
    Ryy += rhs.getRefRyy();
    Rzz += rhs.getRefRzz();
    Rxy += rhs.getRefRxy();
    Rxz += rhs.getRefRxz();
    Ryz += rhs.getRefRyz();
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::timesAssign(ValueType rhs)
{
    VX *= rhs;
    VY *= rhs;
    VZ *= rhs;
    Sxx *= rhs;
    Syy *= rhs;
    Szz *= rhs;
    Sxy *= rhs;
    Sxz *= rhs;
    Syz *= rhs;
    Rxx *= rhs;
    Ryy *= rhs;
    Rzz *= rhs;
    Rxy *= rhs;
    Rxz *= rhs;
    Ryz *= rhs;
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::applyWavefieldTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
{
    VX = lhs * rhs.getRefVX();
    VY = lhs * rhs.getRefVY();
    VZ = lhs * rhs.getRefVZ();
    Sxx = lhs * rhs.getRefSxx();
    Syy = lhs * rhs.getRefSyy();
    Szz = lhs * rhs.getRefSzz();
    Sxy = lhs * rhs.getRefSxy();
    Sxz = lhs * rhs.getRefSxz();
    Syz = lhs * rhs.getRefSyz();
    Rxx = lhs * rhs.getRefRxx();
    Ryy = lhs * rhs.getRefRyy();
    Rzz = lhs * rhs.getRefRzz();
    Rxy = lhs * rhs.getRefRxy();
    Rxz = lhs * rhs.getRefRxz();
    Ryz = lhs * rhs.getRefRyz();
}

template class KITGPI::Wavefields::FD3Dvisco<float>;
template class KITGPI::Wavefields::FD3Dvisco<double>;
