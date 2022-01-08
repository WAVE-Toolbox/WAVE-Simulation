#include "Wavefields3Dviscoelastic.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template <typename ValueType>
scai::hmemo::ContextPtr KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getContextPtr()
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
KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::FD3Dviscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    equationType = "viscoelastic";
    numDimension = 3;
    init(ctx, dist, numRelaxationMechanisms_in);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
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
    
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    scai::lama::DenseVector<ValueType> temp;
    this->initWavefield(temp, ctx, dist);
    Rxx.clear();
    Ryy.clear();
    Rzz.clear();
    Ryz.clear();
    Rxz.clear();
    Rxy.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxx.push_back(temp);
        Ryy.push_back(temp);
        Rzz.push_back(temp);
        Ryz.push_back(temp);
        Rxz.push_back(temp);
        Rxy.push_back(temp);
    }
}

template <typename ValueType>
ValueType KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in)
{
    /* 15 Wavefields in 2D acoustic modeling: Sxx,Syy,Szz,Sxy,Sxz,Syz,Rxx,Ryy,Rzz,Rxy,Rxz,Ryz, Vx, Vy, Vz */
    IndexType numWavefields = 9 + 6 * numRelaxationMechanisms_in;
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
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::write(IndexType snapType, std::string baseName, IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, IndexType fileFormat)
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
        IO::writeVector(Sxx, fileName + ".Sxx." + timeStep, fileFormat);
        IO::writeVector(Syy, fileName + ".Syy." + timeStep, fileFormat);
        IO::writeVector(Szz, fileName + ".Szz." + timeStep, fileFormat);
        IO::writeVector(Sxy, fileName + ".Sxy." + timeStep, fileFormat);
        IO::writeVector(Sxz, fileName + ".Sxz." + timeStep, fileFormat);
        IO::writeVector(Syz, fileName + ".Syz." + timeStep, fileFormat);
        for (int l=0; l<numRelaxationMechanisms; l++) {
            IO::writeVector(Rxx[l], fileName + ".Rxx" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Ryy[l], fileName + ".Ryy" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Rzz[l], fileName + ".Rzz" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Rxy[l], fileName + ".Rxy" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Rxz[l], fileName + ".Rxz" + std::to_string(l+1) + "." + timeStep, fileFormat);
            IO::writeVector(Ryz[l], fileName + ".Ryz" + std::to_string(l+1) + "." + timeStep, fileFormat);
        }
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

/*! \brief decompose wavefields to parts.
 \param decomposeWavefieldType decomposeWavefieldType
 \param wavefieldsDerivative the time derivative of wavefields
 \param derivatives the spatial derivatives
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::decompose(IndexType decomposeWavefieldType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives)
{ 
}

/*! \brief Set all wavefields to zero.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::resetWavefields()
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        this->resetWavefield(Rxx[l]);
        this->resetWavefield(Ryy[l]);
        this->resetWavefield(Rzz[l]);
        this->resetWavefield(Ryz[l]);
        this->resetWavefield(Rxz[l]);
        this->resetWavefield(Rxy[l]);
    }
}

/*! \brief Get numDimension (3)
 */
template <typename ValueType>
int KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getNumDimension() const
{
    return (numDimension);
}

/*! \brief Get equationType (viscoelastic)
 */
template <typename ValueType>
std::string KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getEquationType() const
{
    return (equationType);
}

template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus)
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
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, lama::Vector<ValueType> const &PWaveModulus)
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

//! \brief Not valid in the 3D viscoelastic case
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::getRefP()
{
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D viscoelastic case.")
    return (P);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoelastic<ValueType> KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Wavefields::FD3Dviscoelastic<ValueType> result;
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
KITGPI::Wavefields::FD3Dviscoelastic<ValueType> operator*(ValueType lhs, KITGPI::Wavefields::FD3Dviscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoelastic<ValueType> KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::operator*=(ValueType rhs)
{
    return rhs * *this;
}

/*! \brief Overloading * Operation
 *
 \param rhs seperate Wavefield whith which the components of the current wavefield are multiplied.
 */
template <typename ValueType>
KITGPI::Wavefields::FD3Dviscoelastic<ValueType> KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::operator*(KITGPI::Wavefields::FD3Dviscoelastic<ValueType> rhs)
{
    KITGPI::Wavefields::FD3Dviscoelastic<ValueType> result;
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
KITGPI::Wavefields::FD3Dviscoelastic<ValueType> KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::operator*=(KITGPI::Wavefields::FD3Dviscoelastic<ValueType> rhs)
{
    return rhs * *this;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract wavefield which is assigned.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
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
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxx[l] -= rhs.getRefRxx()[l];
        Ryy[l] -= rhs.getRefRyy()[l];
        Rzz[l] -= rhs.getRefRzz()[l];
        Rxy[l] -= rhs.getRefRxy()[l];
        Rxz[l] -= rhs.getRefRxz()[l];
        Ryz[l] -= rhs.getRefRyz()[l];
    }
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstarct wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs)
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxx[l] += rhs.getRefRxx()[l];
        Ryy[l] += rhs.getRefRyy()[l];
        Rzz[l] += rhs.getRefRzz()[l];
        Rxy[l] += rhs.getRefRxy()[l];
        Rxz[l] += rhs.getRefRxz()[l];
        Ryz[l] += rhs.getRefRyz()[l];
    }
}

/*! \brief function for overloading *= Operation (called in base class)
 *
 \param rhs Scalar which is multiplied.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::timesAssign(ValueType rhs)
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxx[l] *= rhs;
        Ryy[l] *= rhs;
        Rzz[l] *= rhs;
        Rxy[l] *= rhs;
        Rxz[l] *= rhs;
        Ryz[l] *= rhs;
    }
}

/*! \brief apply model transform to wavefields in inversion
 *
 \param rhs Abstract wavefield which is added.
 */
template <typename ValueType>
void KITGPI::Wavefields::FD3Dviscoelastic<ValueType>::applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs)
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Rxx[l] = lhs * rhs.getRefRxx()[l];
        Ryy[l] = lhs * rhs.getRefRyy()[l];
        Rzz[l] = lhs * rhs.getRefRzz()[l];
        Rxy[l] = lhs * rhs.getRefRxy()[l];
        Rxz[l] = lhs * rhs.getRefRxz()[l];
        Ryz[l] = lhs * rhs.getRefRyz()[l];
    }
}

template class KITGPI::Wavefields::FD3Dviscoelastic<float>;
template class KITGPI::Wavefields::FD3Dviscoelastic<double>;
