#include "Viscoelastic.hpp"
using namespace scai;

/*! \brief Switch the default parameterization of this class to modulus
 *
 * This will recalulcate the modulus vectors from the velocity vectors.
 * Moreover, the parametrisation value will be set to zero.
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::switch2modulus()
{
    if (parametrisation == 1) {
        this->calcModuleFromVelocity(velocityP, density, pWaveModulus);
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::switch2velocity()
{
    if (parametrisation == 0) {
        this->calcVelocityFromModule(pWaveModulus, density, velocityP);
        this->calcVelocityFromModule(sWaveModulus, density, velocityS);
        dirtyFlagAveraging = true;
        dirtyFlagModulus = false;
        dirtyFlagVelocity = false;
        parametrisation = 1;
    }
};

/*! \brief Refresh the velocity vectors with the module values
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::refreshVelocity()
{
    if (parametrisation == 0) {
        this->calcVelocityFromModule(pWaveModulus, density, velocityP);
        this->calcVelocityFromModule(sWaveModulus, density, velocityS);
        dirtyFlagVelocity = false;
    }
};

/*! \brief Refresh the module vectors with the velocity values
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::refreshModule()
{
    if (parametrisation == 1) {
        this->calcModuleFromVelocity(velocityP, density, pWaveModulus);
        this->calcModuleFromVelocity(velocityS, density, sWaveModulus);
        dirtyFlagModulus = false;
        dirtyFlagAveraging = true;
    }
};

/*! \brief Prepare modellparameter for visco-elastic modelling
 *
 * Applies Equation 12 from Bohlen 2002 and refreshes the module
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "Preparation of the model parameters…\n");

    refreshModule();

    initializeMatrices(dist, ctx, config, comm);

    /* Set circular frequency w = 2 * pi * relaxation frequency */
    ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
    ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency);

    ValueType sum = w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);

    lama::DenseVector<ValueType> temp(tauS.getDistributionPtr());

    /* Scaling the S-wave Modulus */
    temp = 1.0;
    temp += sum * tauS;
    temp.invert();
    sWaveModulus.scale(temp);

    /* Scaling the P-wave Modulus */
    temp = 1.0;
    temp += sum * tauP;
    temp.invert();
    pWaveModulus.scale(temp);

    calculateAveraging();
    this->getInverseDensity();

    HOST_PRINT(comm, "Model ready!\n\n");
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
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
        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        ValueType getPWaveModulus = config.get<ValueType>("rho") * config.get<ValueType>("velocityP") * config.get<ValueType>("velocityP");
        ValueType getSWaveModulus = config.get<ValueType>("rho") * config.get<ValueType>("velocityS") * config.get<ValueType>("velocityS");
        init(ctx, dist, getPWaveModulus, getSWaveModulus, config.get<ValueType>("rho"), config.get<ValueType>("tauP"), config.get<ValueType>("tauS"), config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));
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
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    init(ctx, dist, pWaveModulus_const, sWaveModulus_const, rho_const, tauP_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in);
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the five given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho_const, scai::lama::Scalar tauP_const, scai::lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    parametrisation = 0;
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    this->initModelparameter(pWaveModulus, ctx, dist, pWaveModulus_const);
    this->initModelparameter(sWaveModulus, ctx, dist, sWaveModulus_const);
    this->initModelparameter(density, ctx, dist, rho_const);
    this->initModelparameter(tauS, ctx, dist, tauS_const);
    this->initModelparameter(tauP, ctx, dist, tauP_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    parametrisation = 0;
    std::string filenamePWaveModulus = filename + ".pWaveModulus.mtx";
    std::string filenameSWaveModulus = filename + ".sWaveModulus.mtx";
    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";

    this->initModelparameter(pWaveModulus, ctx, dist, filenamePWaveModulus, partitionedIn);
    this->initModelparameter(sWaveModulus, ctx, dist, filenameSWaveModulus, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
    this->initModelparameter(tauS, ctx, dist, filenameTauS, partitionedIn);
    this->initModelparameter(tauP, ctx, dist, filenameTauP, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(const Viscoelastic &rhs)
{
    pWaveModulus = rhs.pWaveModulus;
    sWaveModulus = rhs.sWaveModulus;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagModulus = rhs.dirtyFlagModulus;
    dirtyFlagVelocity = rhs.dirtyFlagVelocity;
    parametrisation = rhs.parametrisation;
    inverseDensity = rhs.inverseDensity;
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 *
 *  Calculates pWaveModulus with
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initVelocities(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    parametrisation = 1;
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";

    this->initModelparameter(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initModelparameter(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
    this->initModelparameter(tauS, ctx, dist, filenameTauS, partitionedIn);
    this->initModelparameter(tauP, ctx, dist, filenameTauP, partitionedIn);
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{

    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";

    this->writeModelparameter(density, filenamedensity, partitionedOut);
    this->writeModelparameter(tauP, filenameTauP, partitionedOut);
    this->writeModelparameter(tauS, filenameTauS, partitionedOut);

    std::string filenameP;
    std::string filenameS;

    SCAI_ASSERT_DEBUG(parametrisation == 0 || parametrisation == 1, "Unkown parametrisation");

    switch (parametrisation) {
    case 0:
        filenameP = filename + ".pWaveModulus.mtx";
        filenameS = filename + ".sWaveModulus.mtx";
        this->writeModelparameter(pWaveModulus, filenameP, partitionedOut);
        this->writeModelparameter(sWaveModulus, filenameS, partitionedOut);
    case 1:
        filenameP = filename + ".vp.mtx";
        filenameS = filename + ".vs.mtx";
        this->writeModelparameter(velocityP, filenameP, partitionedOut);
        this->writeModelparameter(velocityS, filenameS, partitionedOut);
        break;
    }
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration config, scai::dmemo::CommunicatorPtr comm)
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, scai::dmemo::CommunicatorPtr /*comm*/)
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calculateAveraging()
{
    this->calculateInverseAveragedDensity(density, inverseDensityAverageX, DensityAverageMatrixX);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageY, DensityAverageMatrixY);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, DensityAverageMatrixZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, sWaveModulusAverageMatrixYZ);
    this->calculateAveragedTauS(tauS, tauSAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedTauS(tauS, tauSAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedTauS(tauS, tauSAverageYZ, sWaveModulusAverageMatrixYZ);
    dirtyFlagAveraging = false;
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-elastic simulation")
    }
    if (relaxationFrequency_in <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-elastic simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*(scai::lama::Scalar rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density * rhs;
    result.tauS = this->tauS * rhs;
    result.tauP = this->tauP * rhs;
    if (parametrisation == 0) {
        result.pWaveModulus = this->pWaveModulus * rhs;
        result.sWaveModulus = this->sWaveModulus * rhs;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityP = this->velocityP * rhs;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> operator*(scai::lama::Scalar lhs, KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*=(scai::lama::Scalar rhs)
{
    return *this * rhs;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density + rhs.density;
    result.tauS = this->tauS + rhs.tauS;
    result.tauP = this->tauP + rhs.tauP;
    if (parametrisation == 0) {
        result.pWaveModulus = this->pWaveModulus + rhs.pWaveModulus;
        result.sWaveModulus = this->sWaveModulus + rhs.sWaveModulus;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityP = this->velocityP + rhs.velocityP;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return *this + rhs;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density - rhs.density;
    result.tauS = this->tauS - rhs.tauS;
    result.tauP = this->tauP - rhs.tauP;
    if (parametrisation == 0) {
        result.pWaveModulus = this->pWaveModulus - rhs.pWaveModulus;
        result.sWaveModulus = this->sWaveModulus - rhs.sWaveModulus;
        return result;
    }
    if (parametrisation == 1) {
        result.velocityP = this->velocityP - rhs.velocityP;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return *this - rhs;
}

template class KITGPI::Modelparameter::Viscoelastic<float>;
template class KITGPI::Modelparameter::Viscoelastic<double>;
