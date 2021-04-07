#include "Elastic.hpp"
#include "../Acquisition/AcquisitionSettings.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 *
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Elastic<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 11 Parameter in elastic modeling: Vp, Vs, rho, invRhoX, invRhoY, invRhoZ, bulk modulus, sWaveModulus, sWaveModulusXY, sWaveModulusXZ, sWaveModulusYZ */
    IndexType numParameter = 11;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modellparameter for modelling
 *
 * Refreshes the modulus, calculates inverseDensity and calculates average values on staggered grid
 *
 \param modelCoordinates coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    // refreshModulus
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
void KITGPI::Modelparameter::Elastic<ValueType>::applyThresholds(Configuration::Configuration const &config)
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

/*! \brief If stream configuration is used, get a subset model from the big model
 \param modelSubset subset model
 \param modelCoordinates coordinate class object of the subset
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinates cut coordinate
 \param cutCoordInd cut coordinate index
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::getModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &modelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd)
{
    auto distBig = velocityP.getDistributionPtr();
    auto dist = modelSubset.getVelocityP().getDistributionPtr();
//     auto comm = dist.getCommunicatorPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist,distBig,modelCoordinates,modelCoordinatesBig,cutCoordinates.at(cutCoordInd));

    lama::DenseVector<ValueType> temp;
    temp = modelSubset.getVelocityP();
    
    temp = shrinkMatrix*velocityP;
    modelSubset.setVelocityP(temp);
        
    temp = shrinkMatrix*velocityS;
    modelSubset.setVelocityS(temp);
//    IO::writeVector(temp, "model/usedSubset_" + std::to_string(cutCoordInd) + ".vs", 2);
    
    temp = shrinkMatrix*density;
    modelSubset.setDensity(temp);
}

/*! \brief If stream configuration is used, set a subset model into the big model
 \param modelSubset subset model
 \param modelCoordinates coordinate class object of the subset
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinates cut coordinate
 \param cutCoordInd cut coordinate index
 \param smoothRange range in x direction which is to be smoothened
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::setModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &invertedModelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NY, scai::IndexType NXBig, scai::IndexType NYBig, scai::IndexType boundaryWidth)
{
    auto distBig = velocityP.getDistributionPtr();
    auto dist = invertedModelSubset.getVelocityP().getDistributionPtr();
//     auto comm = dist.getCommunicatorPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist,distBig,modelCoordinates,modelCoordinatesBig,cutCoordinates.at(cutCoordInd));
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist,distBig,modelCoordinates,modelCoordinatesBig,cutCoordinates.at(cutCoordInd),NX,NYBig,boundaryWidth);
    
    lama::DenseVector<ValueType> temp;
    temp = invertedModelSubset.getVelocityP(); //results of inverted subset model
    //damp the boundary borders
    for (IndexType y = 0; y < NY; y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            temp[modelCoordinates.coordinate2index(i, y, 0)] = temp[modelCoordinates.coordinate2index(i, y, 0)]*(i+1)/boundaryWidth;
            temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)] = temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)]*(i+1)/boundaryWidth;
        }
    }
    temp = shrinkMatrix*temp; //transform subset into big model
    velocityP *= eraseVector;
    velocityP += temp; //take over the values
  
    temp = invertedModelSubset.getVelocityS(); //results of inverted subset model
    //damp the boundary borders
    for (IndexType y = 0; y < NY; y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            temp[modelCoordinates.coordinate2index(i, y, 0)] = temp[modelCoordinates.coordinate2index(i, y, 0)]*(i+1)/boundaryWidth;
            temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)] = temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)]*(i+1)/boundaryWidth;
        }
    }
    temp = shrinkMatrix*temp; //transform subset into big model
    velocityS *= eraseVector;
    velocityS += temp; //take over the values
    
    scai::lama::DenseVector<ValueType> smoothParameter = this->smoothParameter(modelCoordinatesBig, velocityS, cutCoordinates.at(cutCoordInd), smoothRange, NX, NXBig, NYBig);
    velocityS = smoothParameter;
//    IO::writeVector(velocityS, "model/setSubset_" + std::to_string(cutCoordInd) + ".vs" ,2);
    
    temp = invertedModelSubset.getDensity(); //results of inverted subset model
    //damp the boundary borders
    for (IndexType y = 0; y < NY; y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            temp[modelCoordinates.coordinate2index(i, y, 0)] = temp[modelCoordinates.coordinate2index(i, y, 0)]*(i+1)/boundaryWidth;
            temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)] = temp[modelCoordinates.coordinate2index(NX-1-i, y, 0)]*(i+1)/boundaryWidth;
        }
    }
    temp = shrinkMatrix*temp; //transform subset into big model
    density *= eraseVector;
    density += temp; //take over the values

}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "elastic";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2 || config.get<IndexType>("UseVariableGrid") == 1, "Read variable model (ModelRead=2) not available if regular grid is chosen!")

    if ((config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 0) || (config.get<IndexType>("ModelRead") == 2)) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (elastic) parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else if (config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 1) {

        Acquisition::Coordinates<ValueType> regularCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        dmemo::DistributionPtr regularDist(new dmemo::BlockDistribution(regularCoordinates.getNGridpoints(), dist->getCommunicatorPtr()));

        init(ctx, regularDist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "reading regular model finished\n\n")

        init(dist, modelCoordinates, regularCoordinates);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "initialising model on discontineous grid finished\n")

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"));
    }

}

/*! \brief initialisation function which creates a variable grid model on top of a regular model
             \param model regular input model
             \param variableDist Distribution for a variable grid
             \param variableCoordinates Coordinate Class of a Variable Grid
             \param regularCoordinates Coordinate Class of a regular Grid
             */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates)
{

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    variableDist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = variableCoordinates.index2coordinate(ownedIndex);
        IndexType const &regularIndex = regularCoordinates.coordinate2index(coordinate);
        assembly.push(ownedIndex, regularIndex, 1.0);
    }

    lama::CSRSparseMatrix<ValueType> meshingMatrix;
    meshingMatrix = lama::zero<lama::CSRSparseMatrix<ValueType>>(variableDist, velocityP.getDistributionPtr());
    meshingMatrix.fillFromAssembly(assembly);

    density = meshingMatrix * density;
    velocityP = meshingMatrix * velocityP;
    velocityS = meshingMatrix * velocityS;
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
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType pWaveModulus_const, ValueType sWaveModulus_const, ValueType rho)
{
    equationType = "elastic";
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
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho)
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
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "elastic";
    init(ctx, dist, filename, fileFormat);
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    this->initModelparameter(velocityP, ctx, dist, filename + ".vp", fileFormat);
    this->initModelparameter(velocityS, ctx, dist, filename + ".vs", fileFormat);
    this->initModelparameter(density, ctx, dist, filename + ".density", fileFormat);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(const Elastic &rhs)
{
    equationType = rhs.equationType;
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
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::write(std::string filename, IndexType fileFormat) const
{
    IO::writeVector(density, filename + ".density", fileFormat);
    IO::writeVector(velocityP, filename + ".vp", fileFormat);
    IO::writeVector(velocityS, filename + ".vs", fileFormat);
};

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr /*comm*/)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.Elastic.initializeMatrices")

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
void KITGPI::Modelparameter::Elastic<ValueType>::purgeMatrices()
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
void KITGPI::Modelparameter::Elastic<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.Elastic.calculateAveraging")
        this->calculateInverseAveragedDensity(density, inverseDensityAverageX, averageMatrixX);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageY, averageMatrixY);
        this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, averageMatrixZ);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, averageMatrixXY);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, averageMatrixXZ);
        this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, averageMatrixYZ);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Get equationType (elastic)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::Elastic<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to inverse density
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getInverseDensity()
{
    COMMON_THROWEXCEPTION("Inverse density is not set for elastic modelling")
    return (inverseDensity);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return (tauP);
}

/*! \brief Get reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauS() const
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Elastic<ValueType>::getTauSAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return (tauSAverageYZ);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> KITGPI::Modelparameter::Elastic<ValueType>::operator*(ValueType rhs)
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
KITGPI::Modelparameter::Elastic<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::Elastic<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
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
KITGPI::Modelparameter::Elastic<ValueType> &KITGPI::Modelparameter::Elastic<ValueType>::operator=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
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
void KITGPI::Modelparameter::Elastic<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    density = rhs.getDensity();
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
void KITGPI::Modelparameter::Elastic<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP -= rhs.getVelocityP();
    velocityS -= rhs.getVelocityS();
    density -= rhs.getDensity();
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
void KITGPI::Modelparameter::Elastic<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
    density += rhs.getDensity();
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Elastic<float>;
template class KITGPI::Modelparameter::Elastic<double>;
