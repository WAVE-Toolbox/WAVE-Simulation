#include "Elastic.hpp"

using namespace scai;

/*! \brief Prepare modellparameter for modelling
 *
 * Refreshes the modulus, calculates inverseDensity and calculates average values on staggered grid
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "Preparation of the model parameters\n");

    // refreshModulus
    this->getPWaveModulus();
    this->getSWaveModulus();
    initializeMatrices(dist, ctx, config, comm);
    this->getInverseDensity();
    calculateAveraging();

    HOST_PRINT(comm, "Model ready!\n\n");
}



/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::Elastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "Reading model parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"));
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
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho)
{
    init(ctx, dist, pWaveModulus_const, sWaveModulus_const, rho);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar velocityP_const, scai::lama::Scalar velocityS_const, scai::lama::Scalar rho)
{
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the Velocity-Vector "filename".vp.mtx" and "filename".vs.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{  
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initModelparameter(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initModelparameter(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(const Elastic &rhs)
{
    pWaveModulus = rhs.pWaveModulus;
    sWaveModulus = rhs.sWaveModulus;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagPWaveModulus = rhs.dirtyFlagPWaveModulus;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
    inverseDensity = rhs.inverseDensity;
}


/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenameP = filename + ".vp.mtx";
    std::string filenameS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->writeModelparameter(density, filenamedensity, partitionedOut);
    this->writeModelparameter(velocityP, filenameP, partitionedOut);
    this->writeModelparameter(velocityS, filenameS, partitionedOut);
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
void KITGPI::Modelparameter::Elastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration config, scai::dmemo::CommunicatorPtr comm)
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
void KITGPI::Modelparameter::Elastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, scai::dmemo::CommunicatorPtr /*comm*/)
{
    if (dirtyFlagAveraging)
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
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging)
    {
    this->calculateInverseAveragedDensity(density, inverseDensityAverageX, DensityAverageMatrixX);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageY, DensityAverageMatrixY);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, DensityAverageMatrixZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, sWaveModulusAverageMatrixYZ);
    dirtyFlagAveraging = false;
    }
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauS() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Elastic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Elastic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> KITGPI::Modelparameter::Elastic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Modelparameter::Elastic<ValueType> result(*this);
    result *= rhs;
    return result;   
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Modelparameter::Elastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
        density *= rhs;
        velocityP *= rhs;
	velocityS *= rhs;
	
	dirtyFlagInverseDensity=true;
	dirtyFlagPWaveModulus= true;
	dirtyFlagSWaveModulus= true;
	dirtyFlagAveraging= true;
        return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> KITGPI::Modelparameter::Elastic<ValueType>::operator+(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Elastic<ValueType> result(*this);
    result += rhs;
    return result;   
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator+=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    density += rhs.density;
        velocityP += rhs.velocityP;
	velocityS += rhs.velocityS;
	
		dirtyFlagInverseDensity=true;
	dirtyFlagPWaveModulus= true;
	dirtyFlagSWaveModulus= true;
	dirtyFlagAveraging= true;
        return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> KITGPI::Modelparameter::Elastic<ValueType>::operator-(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Elastic<ValueType> result(*this);
    result -= rhs;
    return result; 
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator-=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    density = density -= rhs.density;
        velocityP -= rhs.velocityP;
	velocityS -= rhs.velocityS;
	
		dirtyFlagInverseDensity=true;
	dirtyFlagPWaveModulus= true;
	dirtyFlagSWaveModulus= true;
	dirtyFlagAveraging= true;
        return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    pWaveModulus = rhs.pWaveModulus;
    sWaveModulus = rhs.sWaveModulus;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
	dirtyFlagPWaveModulus= true;
	dirtyFlagSWaveModulus= true;
    inverseDensity = rhs.inverseDensity;
    return *this;
}

template class KITGPI::Modelparameter::Elastic<float>;
template class KITGPI::Modelparameter::Elastic<double>;
