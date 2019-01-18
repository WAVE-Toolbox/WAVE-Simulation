#include "SH.hpp"

using namespace scai;

/*! \brief Prepare modellparameter for modelling
 *
 * Refreshes the modulus, calculates inverse density and average Values on staggered grid
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters…\n");

    // refreshModulus();
    this->getSWaveModulus();
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    this->getInverseDensity();
    calculateAveraging();
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    lama::DenseVector<ValueType> mask(velocityS); //mask to restore vacuum
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);

    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("lowerVPTh"), config.get<ValueType>("lowerVPTh"), 1);
    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("upperVPTh"), config.get<ValueType>("upperVPTh"), 2);
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty

    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("lowerDensityTh"), config.get<ValueType>("lowerDensityTh"), 1);
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("upperDensityTh"), config.get<ValueType>("upperDensityTh"), 2);
    dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed
    dirtyFlagAveraging = true;      // If S-Wave velocity will be changed, averaging needs to be redone

    velocityS *= mask;
    density *= mask;
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "sh";
    init(config, ctx, dist);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model parameter (SH) from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityS"), config.get<ValueType>("rho"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given values.
 \param ctx Context
 \param dist Distribution
 \param sWaveModulus_const S-wave modulus given as 
 \param rho Density given as ValueType
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho)
{
    equationType = "sh";
    init(ctx, dist, velocityS_const, rho);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given values.
 \param ctx Context
 \param dist Distribution
 \param sWaveModulus_const S-wave modulus given as ValueType
 \param rho Density given as ValueType
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const)
{
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the S-wave modulus  ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    equationType = "sh";
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the S-wave velocity ".vs.mtx" and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initModelparameter(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(const SH &rhs)
{
    equationType = rhs.equationType;
    sWaveModulus = rhs.sWaveModulus;
    velocityS = rhs.velocityS;
    density = rhs.density;
    inverseDensity = rhs.inverseDensity;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
}

/*! \brief Write model to an external file
 *
 \param filename For the S-wave velocity ".vs.mtx" and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenameS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->writeModelparameter(density, filenamedensity, partitionedOut);
    this->writeModelparameter(velocityS, filenameS, partitionedOut);
};

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr /*comm*/)
{

    SCAI_REGION("initializeMatrices")
    // reuse of density average for the swavemodulus : sxz and syz are on the same spot as vx and vy in acoustic modeling
    this->calcDensityAverageMatrixX(modelCoordinates, dist);
    this->calcDensityAverageMatrixY(modelCoordinates, dist);

    sWaveModulusAverageMatrixXZ.setContextPtr(ctx);
    sWaveModulusAverageMatrixYZ.setContextPtr(ctx);
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        // reuse of density average for the swavemodulus : sxz and syz are on the same spot as vx and vy in acoustic modeling
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, DensityAverageMatrixX);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, DensityAverageMatrixY);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Get equationType (sh)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::SH<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getVelocityP() const
{
    COMMON_THROWEXCEPTION("There is no velocityP parameter in an sh modelling")
    return (velocityP);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getPWaveModulus() const
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus parameter in an sh modelling")
    return (pWaveModulus);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getTauS() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::SH<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an sh modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::SH<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an sh modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageYZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs ValueType factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs ValueType factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs valueType factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> &KITGPI::Modelparameter::SH<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    velocityS *= rhs;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator+(KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> &KITGPI::Modelparameter::SH<ValueType>::operator+=(KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    density += rhs.density;
    velocityS += rhs.velocityS;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator-(KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> &KITGPI::Modelparameter::SH<ValueType>::operator-=(KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    density = density -= rhs.density;
    velocityS -= rhs.velocityS;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> &KITGPI::Modelparameter::SH<ValueType>::operator=(KITGPI::Modelparameter::SH<ValueType> const &rhs)
{
    sWaveModulus = rhs.sWaveModulus;
    velocityS = rhs.velocityS;
    density = rhs.density;
    inverseDensity = rhs.inverseDensity;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
    dirtyFlagAveraging = rhs.dirtyFlagAveraging;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    sWaveModulus = rhs.getSWaveModulus();
    velocityS = rhs.getVelocityS();
    inverseDensity = rhs.getInverseDensity();
    density = rhs.getDensity();
    dirtyFlagInverseDensity = rhs.getDirtyFlagInverseDensity();
    dirtyFlagSWaveModulus = rhs.getDirtyFlagSWaveModulus();
    dirtyFlagAveraging = rhs.getDirtyFlagAveraging();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityS -= rhs.getVelocityS();
    density -= rhs.getDensity();
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityS += rhs.getVelocityS();
    density += rhs.getDensity();
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}
template class KITGPI::Modelparameter::SH<float>;
template class KITGPI::Modelparameter::SH<double>;
