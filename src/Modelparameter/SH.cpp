#include "SH.hpp"

using namespace scai;

/*! \brief Prepare modellparameter for modelling
 *
 * Refreshes the module if parameterisation is in terms of velocities
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "Preparation of the model parameters…\n");

    refreshModule();
    initializeMatrices(dist, ctx, config, comm);
    this->getInverseDensity();
    calculateAveraging();

    HOST_PRINT(comm, "Model ready!\n\n");
}

/*! \brief Switch the default parameterization of this class to modulus
 *
 * This will recalulcate the modulus vectors from the velocity vectors.
 * Moreover, the parametrisation value will be set to zero.
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::switch2modulus()
{
    if (parametrisation == 1) {
        this->calcModuleFromVelocity(velocityS, density, sWaveModulus);
        dirtyFlagAveraging = true;
        dirtyFlagModulus = false;
        dirtyFlagVelocity = false;
        parametrisation = 0;
    }
};

/*! \brief Switch the default parameterization of this class to velocity
 *
 * This will recalulcate the velocity vectors from the modulus vectors.
 * Moreover, the parametrisation value will be set to one.
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::switch2velocity()
{
    if (parametrisation == 0) {
        this->calcVelocityFromModule(sWaveModulus, density, velocityS);
        dirtyFlagModulus = false;
        dirtyFlagVelocity = false;
        parametrisation = 1;
    }
};

/*! \brief Refresh the velocity vectors with the module values
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::refreshVelocity()
{
    if (parametrisation == 0) {
        this->calcVelocityFromModule(sWaveModulus, density, velocityS);
        dirtyFlagVelocity = false;
    }
};

/*! \brief Refresh the module vectors with the velocity values
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::refreshModule()
{
    if (parametrisation == 1) {
        this->calcModuleFromVelocity(velocityS, density, sWaveModulus);
        dirtyFlagModulus = false;
        dirtyFlagAveraging = true;
    }
};

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
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

        HOST_PRINT(dist->getCommunicatorPtr(), "Reading model parameter from file...\n");

        switch (config.get<IndexType>("ModelParametrisation")) {
        case 1:
            init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));
            break;
        case 2:
            initVelocities(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));
            break;
        default:
            COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
            break;
        }

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        ValueType getSWaveModulus = config.get<ValueType>("rho") * config.get<ValueType>("velocityS") * config.get<ValueType>("velocityS");
        init(ctx, dist, getSWaveModulus, config.get<ValueType>("rho"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho)
{
    init(ctx, dist, sWaveModulus_const, rho);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho)
{
    parametrisation = 0;
    this->initModelparameter(sWaveModulus, ctx, dist, sWaveModulus_const);
    this->initModelparameter(density, ctx, dist, rho);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenameSWaveModulus Name of file that will be read for the S-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filenameSWaveModulus, std::string filenamerho, IndexType partitionedIn)
{
    init(ctx, dist, filenameSWaveModulus, filenamerho, partitionedIn);
}

/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenameSWaveModulus Name of file that will be read for the S-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filenameSWaveModulus, std::string filenamerho, IndexType partitionedIn)
{
    parametrisation = 0;
    this->initModelparameter(density, ctx, dist, filenamerho, partitionedIn);
    this->initModelparameter(sWaveModulus, ctx, dist, filenameSWaveModulus, partitionedIn);
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
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the S-wave modulus ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    parametrisation = 0;
    std::string filenameSWaveModulus = filename + ".sWaveModulus.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initModelparameter(sWaveModulus, ctx, dist, filenameSWaveModulus, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType>::SH(const SH &rhs)
{
    sWaveModulus = rhs.sWaveModulus;
    velocityS = rhs.velocityS;
    density = rhs.density;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagModulus = rhs.dirtyFlagModulus;
    dirtyFlagVelocity = rhs.dirtyFlagVelocity;
    parametrisation = rhs.parametrisation;
    inverseDensity = rhs.inverseDensity;
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the Velocity-Vector "filename".vs.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::initVelocities(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    parametrisation = 1;
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initModelparameter(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
}

/*! \brief Write model to an external file
 *
 \param filenameS Filename for S-wave modulus / S-wave velocity model
 \param filenamedensity Filename for Density model
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::write(std::string filenameS, std::string filenamedensity, IndexType partitionedOut) const
{
    SCAI_ASSERT_DEBUG(parametrisation == 0 || parametrisation == 1, "Unkown parametrisation");

    this->writeModelparameter(density, filenamedensity, partitionedOut);

    switch (parametrisation) {
    case 0:
        this->writeModelparameter(sWaveModulus, filenameS, partitionedOut);
        break;
    case 1:
        this->writeModelparameter(velocityS, filenameS, partitionedOut);
        break;
    }
};

/*! \brief Write model to an external file
 *
 \param filename For the S-wave modulus ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    SCAI_ASSERT_DEBUG(parametrisation == 0 || parametrisation == 1, "Unkown parametrisation");

    std::string filenameS;
    std::string filenamedensity = filename + ".density.mtx";

    switch (parametrisation) {
    case 0:
        filenameS = filename + ".sWaveModulus.mtx";
        break;
    case 1:
        filenameS = filename + ".vs.mtx";
        break;
    }

    write(filenameS, filenamedensity, partitionedOut);
};

//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration config, scai::dmemo::CommunicatorPtr comm)
{
    initializeMatrices(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"), config.get<ValueType>("DT"), comm);
}

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, scai::dmemo::CommunicatorPtr /*comm*/)
{

    SCAI_REGION("initializeMatrices")

    this->calcDensityAverageMatrixX(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixY(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixZ(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixXY(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixXZ(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixYZ(NX, NY, NZ, dist);

    DensityAverageMatrixX.setContextPtr(ctx);
    DensityAverageMatrixY.setContextPtr(ctx);
    DensityAverageMatrixZ.setContextPtr(ctx);
    sWaveModulusAverageMatrixXY.setContextPtr(ctx);
    sWaveModulusAverageMatrixXZ.setContextPtr(ctx);
    sWaveModulusAverageMatrixYZ.setContextPtr(ctx);
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::SH<ValueType>::calculateAveraging()
{

    this->calculateInverseAveragedDensity(density, inverseDensityAverageX, DensityAverageMatrixX);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageY, DensityAverageMatrixY);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, DensityAverageMatrixZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, sWaveModulusAverageMatrixYZ);
    dirtyFlagAveraging = false;
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getVelocityP()
{
    COMMON_THROWEXCEPTION("There is no velocityP parameter in an sh modelling")
    return (velocityP);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getPWaveModulus()
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus parameter in an sh modelling")
    return (pWaveModulus);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getTauP()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getTauS()
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
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::SH<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageYZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result;
    result.density = this->density * rhs;
    if (parametrisation == 0) {
        result.sWaveModulus = this->sWaveModulus * rhs;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityS = this->velocityS * rhs;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Modelparameter::SH<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return *this * rhs;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator+(KITGPI::Modelparameter::SH<ValueType> rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result;
    result.density = this->density + rhs.density;
    if (parametrisation == 0) {
        result.sWaveModulus = this->sWaveModulus + rhs.sWaveModulus;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityS = this->velocityS + rhs.velocityS;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator+=(KITGPI::Modelparameter::SH<ValueType> rhs)
{
    return *this + rhs;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator-(KITGPI::Modelparameter::SH<ValueType> rhs)
{
    KITGPI::Modelparameter::SH<ValueType> result;
    result.density = this->density - rhs.density;
    if (parametrisation == 0) {
        result.sWaveModulus = this->sWaveModulus - rhs.sWaveModulus;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityS = this->velocityS - rhs.velocityS;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::SH<ValueType> KITGPI::Modelparameter::SH<ValueType>::operator-=(KITGPI::Modelparameter::SH<ValueType> rhs)
{
    return *this - rhs;
}

template class KITGPI::Modelparameter::SH<float>;
template class KITGPI::Modelparameter::SH<double>;
