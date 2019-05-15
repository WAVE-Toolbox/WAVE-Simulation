#include "Acoustic.hpp"
#include <scai/lama/io/FileIO.hpp>
using namespace scai;
using namespace KITGPI;

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
void KITGPI::Modelparameter::Acoustic<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    // refreshModulus();
    this->getPWaveModulus();
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    this->getInverseDensity();
    calculateAveraging();
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    lama::DenseVector<ValueType> mask(velocityP); //mask to restore vacuum
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);

    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("lowerVPTh"), config.get<ValueType>("lowerVPTh"), 1);
    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("upperVPTh"), config.get<ValueType>("upperVPTh"), 2);
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty

    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("lowerDensityTh"), config.get<ValueType>("lowerDensityTh"), 1);
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("upperDensityTh"), config.get<ValueType>("upperDensityTh"), 2);
    dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed
    dirtyFlagAveraging = true;      // If S-Wave velocity will be changed, averaging needs to be redone

    velocityP *= mask;
    density *= mask;
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "acoustic";
    init(config, ctx, dist);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model parameter (acoustic) from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("rho"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
    }
}

/*! \brief initialisation function which creates a variable grid model on top of a regular model
             \param config Configuration from configuration file
             \param ctx Context
             \param variableDist Distribution for a variable grid
             \param regularDist Distribution for a regular grid
             \param variableCoordinates Coordinate Class of a Variable Grid
             \param regularCoordinates Coordinate Class of a regular Grid
             */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(KITGPI::Modelparameter::Modelparameter<ValueType> const &model, scai::dmemo::DistributionPtr variableDist,Acquisition::Coordinates<ValueType> const &variableCoordinates,Acquisition::Coordinates<ValueType> const &regularCoordinates)
{
   auto const &regularDensity = model.getDensity();
   auto const &regularVelocityP = model.getVelocityP();
   this->initModelparameter(density, regularDensity.getContextPtr(), variableDist, 0);
   this->initModelparameter(velocityP, regularDensity.getContextPtr(), variableDist, 0);
   
   
    for (IndexType variableIndex=0;variableIndex<variableCoordinates.getNGridpoints();variableIndex++) {
        Acquisition::coordinate3D coordinate = variableCoordinates.index2coordinate(variableIndex);
        IndexType const &regularIndex=regularCoordinates.coordinate2index(coordinate);
        density.setValue(variableIndex,regularDensity.getValue(regularIndex));
        velocityP.setValue(variableIndex,regularVelocityP.getValue(regularIndex));
    }

}
/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType pWaveModulus_const, ValueType rho_const)
{
    equationType = "acoustic";
    init(ctx, dist, pWaveModulus_const, rho_const);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType rho_const)
{
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(density, ctx, dist, rho_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    equationType = "acoustic";
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Velocity-Vector "filename".vp.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->initModelparameter(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(const Acoustic &rhs)
{
    equationType = rhs.equationType;
    pWaveModulus = rhs.pWaveModulus;
    velocityP = rhs.velocityP;
    inverseDensity = rhs.inverseDensity;
    density = rhs.density;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagPWaveModulus = rhs.dirtyFlagPWaveModulus;
}

/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenameP = filename + ".vp.mtx";
    std::string filenamedensity = filename + ".density.mtx";

    this->writeModelparameter(density, filenamedensity, partitionedOut);
    this->writeModelparameter(velocityP, filenameP, partitionedOut);
};

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr /*comm*/)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("initializeMatrices")

        this->calcDensityAverageMatrixX(modelCoordinates, dist);
        this->calcDensityAverageMatrixY(modelCoordinates, dist);
        this->calcDensityAverageMatrixZ(modelCoordinates, dist);

        DensityAverageMatrixX.setContextPtr(ctx);
        DensityAverageMatrixY.setContextPtr(ctx);
        DensityAverageMatrixZ.setContextPtr(ctx);
    }
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        this->calculateInverseAveragedDensity(density, inverseDensityAverageX, DensityAverageMatrixX);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageY, DensityAverageMatrixY);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, DensityAverageMatrixZ);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Get equationType (acoustic)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::Acoustic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to S-wave modulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulus()
{
    COMMON_THROWEXCEPTION("S-wave modulus is not set for acoustic modelling")
    return (sWaveModulus);
}

/*! \brief Get reference to S-wave modulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulus() const
{
    COMMON_THROWEXCEPTION("S-wave modulus is not set for acoustic modelling")
    return (sWaveModulus);
}

/*! \brief Get reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getVelocityS() const
{
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return (velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauS() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauS);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Acoustic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Acoustic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Get reference to S-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXY()
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageXY);
}

/*! \brief Get reference to S-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXY() const
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageXY);
}

/*! \brief Get reference to S-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXZ()
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageXZ);
}

/*! \brief Get reference to S-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXZ() const
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageXZ);
}

/*! \brief Get reference to S-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageYZ()
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageYZ);
}

/*! \brief Get reference to S-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageYZ() const
{
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return (sWaveModulusAverageYZ);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief non-member function to multiply (scalar as left operand)
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> &KITGPI::Modelparameter::Acoustic<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    velocityP *= rhs;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator+(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> &KITGPI::Modelparameter::Acoustic<ValueType>::operator+=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs)
{
    density += rhs.density;
    velocityP += rhs.velocityP;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator-(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> &KITGPI::Modelparameter::Acoustic<ValueType>::operator-=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityP -= rhs.velocityP;

    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> &KITGPI::Modelparameter::Acoustic<ValueType>::operator=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs)
{
    pWaveModulus = rhs.pWaveModulus;
    velocityP = rhs.velocityP;
    inverseDensity = rhs.inverseDensity;
    density = rhs.density;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagPWaveModulus = rhs.dirtyFlagPWaveModulus;
    dirtyFlagAveraging = rhs.dirtyFlagAveraging;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    pWaveModulus = rhs.getPWaveModulus();
    velocityP = rhs.getVelocityP();
    inverseDensity = rhs.getInverseDensity();
    density = rhs.getDensity();
    dirtyFlagInverseDensity = rhs.getDirtyFlagInverseDensity();
    dirtyFlagPWaveModulus = rhs.getDirtyFlagPWaveModulus();
    dirtyFlagAveraging = rhs.getDirtyFlagAveraging();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP -= rhs.getVelocityP();
    density -= rhs.getDensity();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP += rhs.getVelocityP();
    density += rhs.getDensity();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Acoustic<double>;
template class KITGPI::Modelparameter::Acoustic<float>;
