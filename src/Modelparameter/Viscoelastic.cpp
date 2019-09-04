#include "Viscoelastic.hpp"
#include "../IO/IO.hpp"
using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 * 
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Viscoelastic<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 15 Parameter in Viscoelastic modeling:  rho, Vp, Vs, invRhoX,invRhoY, invRhoZ,  bulk modulus, sWaveModulus, sWaveModulusXY, sWaveModulusXZ, sWaveModulus YZ, tauP, tauS, tauSXY, tauSXZ, tauSYZ*/
    IndexType numParameter = 15;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modellparameter for visco-elastic modelling
 *
 * Applies Equation 12 from Bohlen 2002 and refreshes the modulus
 *
 \param modelCoordinates coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    //refreshModulus
    this->getPWaveModulus();
    this->getSWaveModulus();
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    calculateAveraging();
    purgeMatrices();
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    lama::DenseVector<ValueType> maskP(velocityP); //mask to restore vacuum
    maskP.unaryOp(maskP, common::UnaryOp::SIGN);
    maskP.unaryOp(maskP, common::UnaryOp::ABS);

    lama::DenseVector<ValueType> maskS(velocityS); //mask to restore acoustic media
    maskS.unaryOp(maskP, common::UnaryOp::SIGN);
    maskS.unaryOp(maskP, common::UnaryOp::ABS);

    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("lowerVPTh"), config.get<ValueType>("lowerVPTh"), 1);
    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("upperVPTh"), config.get<ValueType>("upperVPTh"), 2);
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty

    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("lowerDensityTh"), config.get<ValueType>("lowerDensityTh"), 1);
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("upperDensityTh"), config.get<ValueType>("upperDensityTh"), 2);
    dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed

    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("lowerVSTh"), config.get<ValueType>("lowerVSTh"), 1);
    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("upperVSTh"), config.get<ValueType>("upperVSTh"), 2);
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone

    velocityP *= maskP;
    density *= maskP;
    velocityS *= maskS;
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "viscoelastic";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    if (config.get<IndexType>("ModelRead") == 1) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (viscoelastic) parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));
        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else if (config.get<IndexType>("ModelRead") == 2) {

        Acquisition::Coordinates<ValueType> regularCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        dmemo::DistributionPtr regularDist(new dmemo::BlockDistribution(regularCoordinates.getNGridpoints(), dist->getCommunicatorPtr()));

        init(ctx, regularDist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "reading regular model finished\n\n")

        init(dist, modelCoordinates, regularCoordinates);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "initialising model on discontineous grid finished\n")

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"), config.get<ValueType>("tauP"), config.get<ValueType>("tauS"), config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("FileFormat"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityP_const P-wave velocity given as Scalar
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    equationType = "viscoelastic";
    init(ctx, dist, velocityP_const, velocityS_const, rho_const, tauP_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in);
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
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
 \param fileFormat Input file format 0=mtx 1=lmf
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "viscoelastic";
    init(ctx, dist, filename, fileFormat);
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param fileFormat Input file format 0=mtx 1=lmf
 *
 *  Calculates pWaveModulus with
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    this->initModelparameter(velocityS, ctx, dist, filename + ".vs", fileFormat);
    this->initModelparameter(velocityP, ctx, dist, filename + ".vp", fileFormat);
    this->initModelparameter(density, ctx, dist, filename + ".density", fileFormat);
    this->initModelparameter(tauS, ctx, dist, filename + ".tauS", fileFormat);
    this->initModelparameter(tauP, ctx, dist, filename + ".tauP", fileFormat);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(const Viscoelastic &rhs)
{
    equationType = rhs.equationType;
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
    dirtyFlagPWaveModulus = rhs.dirtyFlagPWaveModulus;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
    inverseDensity = rhs.inverseDensity;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave velocity ".vp.mtx" is added, for the S-wave velocity ".vs.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param fileFormat Output file format 0=mtx 1=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::write(std::string filename, IndexType fileFormat) const
{
    IO::writeVector(density, filename + ".density", fileFormat);
    IO::writeVector(tauP, filename + ".tauP", fileFormat);
    IO::writeVector(tauS, filename + ".tauS", fileFormat);
    IO::writeVector(velocityP, filename + ".vp", fileFormat);
    IO::writeVector(velocityS, filename + ".vs", fileFormat);
};

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates coordinate class object
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr /*comm*/)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("initializeMatrices")

        this->calcAverageMatrixX(modelCoordinates, dist);
        this->calcAverageMatrixY(modelCoordinates, dist);
        this->calcAverageMatrixZ(modelCoordinates, dist);
        this->calcAverageMatrixXY(modelCoordinates, dist);
        this->calcAverageMatrixXZ(modelCoordinates, dist);
        this->calcAverageMatrixYZ(modelCoordinates, dist);

        averageMatrixX.setContextPtr(ctx);
        averageMatrixY.setContextPtr(ctx);
        averageMatrixZ.setContextPtr(ctx);
        averageMatrixXY.setContextPtr(ctx);
        averageMatrixXZ.setContextPtr(ctx);
        averageMatrixYZ.setContextPtr(ctx);
    }
}

//! \brief Purge Averaging matrices to free memory
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::purgeMatrices()
{
    averageMatrixX.purge();
    averageMatrixY.purge();
    averageMatrixZ.purge();
    averageMatrixXY.purge();
    averageMatrixXZ.purge();
    averageMatrixYZ.purge();
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        this->calculateInverseAveragedDensity(density, inverseDensityAverageX, averageMatrixX);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageY, averageMatrixY);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, averageMatrixZ);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, averageMatrixXY);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, averageMatrixXZ);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, averageMatrixYZ);
        this->calculateAveragedTauS(tauS, tauSAverageXY, averageMatrixXY);
        this->calculateAveragedTauS(tauS, tauSAverageXZ, averageMatrixXZ);
        this->calculateAveragedTauS(tauS, tauSAverageYZ, averageMatrixYZ);
        dirtyFlagAveraging = false;
    }
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

/*! \brief Get equationType (viscoelastic)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::Viscoelastic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to inverse density
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Viscoelastic<ValueType>::getInverseDensity()
{
    COMMON_THROWEXCEPTION("Inverse density is not set for elastic modelling")
    return (inverseDensity);
}

/*! \brief Get const reference to P-wave modulus (viscoelastic case)
 *
 * if P-Wave Modulus is dirty (eg. because of changes in velocityP, the P-Wave modulus will be recalculated
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Viscoelastic<ValueType>::getPWaveModulus()
{
    // If the modulus is dirty, than recalculate
    if (dirtyFlagPWaveModulus) {
        HOST_PRINT(velocityP.getDistributionPtr()->getCommunicatorPtr(), "P-Wave modulus will be calculated from density,velocityP,tauP and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityP, density, pWaveModulus);
        /* Set circular frequency w = 2 * pi * relaxation frequency */
        ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
        ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency);

        ValueType sum = w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);

        /* Scaling the P-wave Modulus */

        auto temp = lama::eval<lama::DenseVector<ValueType>>(1.0 + sum * tauP);
        pWaveModulus = pWaveModulus / temp;
        dirtyFlagPWaveModulus = false;
    }

    return (pWaveModulus);
}

/*! \brief Get const reference to P-wave modulus (viscoelastic case)
 *
 * if P-Wave Modulus is dirty (eg. because of changes in velocityP, the P-Wave modulus will be recalculated
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Viscoelastic<ValueType>::getPWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagPWaveModulus == false, "P-Wave Modulus has to be recalculated! ");
    return (pWaveModulus);
}

/*! \brief Get const reference to S-wave modulus (viscoelastic case)
 *
 * if S-Wave Modulus is dirty (eg. because of changes in velocityS, the S-Wave modulus will be recalculated
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Viscoelastic<ValueType>::getSWaveModulus()
{
    // If the modulus is dirty, than recalculate
    if (dirtyFlagSWaveModulus) {
        HOST_PRINT(velocityS.getDistributionPtr()->getCommunicatorPtr(), "S-Wave modulus will be calculated from density,velocityS,tauS and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityS, density, sWaveModulus);
        /* Set circular frequency w = 2 * pi * relaxation frequency */
        ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
        ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency);

        ValueType sum = w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);

        /* Scaling the S-wave Modulus */
        auto temp = lama::eval<lama::DenseVector<ValueType>>(1.0 + sum * tauS);
        sWaveModulus = sWaveModulus / temp;
        dirtyFlagSWaveModulus = false;
    }

    return (sWaveModulus);
}

/*! \brief Get const reference to S-wave modulus (viscoelastic case) const
 *
 * if S-Wave Modulus is dirty (eg. because of changes in velocityS it ca
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Viscoelastic<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "P-Wave Modulus has to be recalculated! ");
    return (sWaveModulus);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    tauS *= rhs;
    tauP *= rhs;
    velocityP *= rhs;
    velocityS *= rhs;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+=(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    density += rhs.density;
    tauS += rhs.tauS;
    tauP += rhs.tauP;
    velocityP += rhs.velocityP;
    velocityS += rhs.velocityS;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-=(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    density = density -= rhs.density;
    tauS -= rhs.tauS;
    tauP -= rhs.tauP;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator=(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    density = rhs.getDensity();
    tauS = rhs.getTauS();
    tauP = rhs.getTauP();
    relaxationFrequency = rhs.getRelaxationFrequency();
    numRelaxationMechanisms = rhs.getNumRelaxationMechanisms();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP -= rhs.getVelocityP();
    velocityS -= rhs.getVelocityS();
    density -= rhs.getDensity();
    tauS -= rhs.getTauS();
    tauP -= rhs.getTauP();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
    density += rhs.getDensity();
    tauS += rhs.getTauS();
    tauP += rhs.getTauP();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Viscoelastic<float>;
template class KITGPI::Modelparameter::Viscoelastic<double>;
