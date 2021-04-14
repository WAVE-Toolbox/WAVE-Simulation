#include "EMEM.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 * 
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::EMEM<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 9 Parameter in emem modeling: inverseMagneticPermeabilityEMXY, inverseMagneticPermeabilityEMYZ, inverseMagneticPermeabilityEMXY, CaAverageX, CaAverageY, CaAverageZ, CbAverageX, CbAverageY, CbAverageZ */
    IndexType numParameter = 9;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modelparameter for modelling
 *
 * Refreshes the modulus, calculates and calculates average values on staggered grid
 *
 \param modelCoordinates coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    // refreshModulus
    this->getVelocityEM();
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    calculateAveraging();
    purgeMatrices();
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    dielectricPermittivityEM /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivityEM
    
    lama::DenseVector<ValueType> mask; //mask to restore vacuum
    mask = dielectricPermittivityEM-1;
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);
    
    Common::searchAndReplace<ValueType>(conductivityEM, config.get<ValueType>("lowerSigmaEMTh"), config.get<ValueType>("lowerSigmaEMTh"), 1);
    Common::searchAndReplace<ValueType>(conductivityEM, config.get<ValueType>("upperSigmaEMTh"), config.get<ValueType>("upperSigmaEMTh"), 2);    
    Common::searchAndReplace<ValueType>(dielectricPermittivityEM, config.get<ValueType>("lowerEpsilonEMrTh"), config.get<ValueType>("lowerEpsilonEMrTh"), 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityEM, config.get<ValueType>("upperEpsilonEMrTh"), config.get<ValueType>("upperEpsilonEMrTh"), 2);
       
    dirtyFlagVelocivityEM = true ;   // If EM-parameters will be changed, velocityEM needs to be redone
    dirtyFlagAveraging = true;      // If EM-parameters will be changed, averaging needs to be redone
      
    Common::searchAndReplace<ValueType>(porosity, config.get<ValueType>("lowerPorosityTh"), config.get<ValueType>("lowerPorosityTh"), 1);
    Common::searchAndReplace<ValueType>(porosity, config.get<ValueType>("upperPorosityTh"), config.get<ValueType>("upperPorosityTh"), 2);
    Common::searchAndReplace<ValueType>(saturation, config.get<ValueType>("lowerSaturationTh"), config.get<ValueType>("lowerSaturationTh"), 1);
    Common::searchAndReplace<ValueType>(saturation, config.get<ValueType>("upperSaturationTh"), config.get<ValueType>("upperSaturationTh"), 2);
    
    conductivityEM *= mask;
    porosity *= mask;
    saturation *= mask;
    
    dielectricPermittivityEM -= 1;
    dielectricPermittivityEM *= mask;
    dielectricPermittivityEM += 1;
    dielectricPermittivityEM *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM   
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivityEM().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);

    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * dielectricPermittivityEM;
    modelPerShot.setDielectricPermittivityEM(temp);
        
    temp = shrinkMatrix * conductivityEM;
    modelPerShot.setConductivityEM(temp);
    
    temp = shrinkMatrix * porosity;
    modelPerShot.setPorosity(temp);
    
    temp = shrinkMatrix * saturation;
    modelPerShot.setSaturation(temp);
}

/*! \brief If stream configuration is used, set a pershot model into the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivityEM().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate, boundaryWidth);
    scai::lama::SparseVector<ValueType> restoreVector;
    restoreVector = 1.0 - eraseVector;
    
    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * modelPerShot.getDielectricPermittivityEM(); //transform pershot into big model
    temp *= restoreVector;
    dielectricPermittivityEM *= eraseVector;
    dielectricPermittivityEM += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getConductivityEM(); //transform pershot into big model
    temp *= restoreVector;
    conductivityEM *= eraseVector;
    conductivityEM += temp; //take over the values
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "emem";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{    
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2 || config.get<IndexType>("UseVariableGrid") == 1, "Read variable model (ModelRead=2) not available if regular grid is chosen!")

    if ((config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 0) || (config.get<IndexType>("ModelRead") == 2)) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (emem) parameter from file...\n");

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
        init(ctx, dist, config.get<ValueType>("muEMr"), config.get<ValueType>("sigmaEM"), config.get<ValueType>("epsilonEMr"));
    }

}

/*! \brief initialisation function which creates a variable grid model on top of a regular model
             \param model regular input model
             \param variableDist Distribution for a variable grid
             \param variableCoordinates Coordinate Class of a Variable Grid
             \param regularCoordinates Coordinate Class of a regular Grid
             */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates)
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
    meshingMatrix = lama::zero<lama::CSRSparseMatrix<ValueType>>(variableDist, conductivityEM.getDistributionPtr());
    meshingMatrix.fillFromAssembly(assembly);

    magneticPermeabilityEM = meshingMatrix * magneticPermeabilityEM;
    conductivityEM = meshingMatrix * conductivityEM;
    dielectricPermittivityEM = meshingMatrix * dielectricPermittivityEM;
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeabilityEM_const magneticPermeabilityEM given as Scalar
 \param conductivityEM_const conductivityEM given as Scalar
 \param dielectricPermittivityEM_const dielectricPermittivityEM given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const)
{
    equationType = "emem";
    init(ctx, dist, magneticPermeabilityEM_const, conductivityEM_const, dielectricPermittivityEM_const);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeabilityEM_const magneticPermeabilityEM given as Scalar
 \param conductivityEM_const conductivityEM given as Scalar
 \param dielectricPermittivityEM_const dielectricPermittivityEM given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const)
{
    magneticPermeabilityEM_const *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeabilityEM
    dielectricPermittivityEM_const *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM
    
    this->initModelparameter(magneticPermeabilityEM, ctx, dist, magneticPermeabilityEM_const);
    this->initModelparameter(conductivityEM, ctx, dist, conductivityEM_const);
    this->initModelparameter(dielectricPermittivityEM, ctx, dist, dielectricPermittivityEM_const);
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
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "emem";
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
void KITGPI::Modelparameter::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    this->initModelparameter(magneticPermeabilityEM, ctx, dist, filename + ".muEMr", fileFormat);
    this->initModelparameter(conductivityEM, ctx, dist, filename + ".sigmaEM", fileFormat);
    this->initModelparameter(dielectricPermittivityEM, ctx, dist, filename + ".epsilonEMr", fileFormat);
    this->initModelparameter(porosity, ctx, dist, filename + ".porosity", fileFormat);
    this->initModelparameter(saturation, ctx, dist, filename + ".saturation", fileFormat);
    
    magneticPermeabilityEM *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeabilityEM
    dielectricPermittivityEM *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(const EMEM &rhs)
{
    equationType = rhs.equationType;
    velocivityEM = rhs.velocivityEM;
    magneticPermeabilityEM = rhs.magneticPermeabilityEM;
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    dirtyFlagVelocivityEM = rhs.dirtyFlagVelocivityEM;
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Write model to an external file
 *
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    scai::lama::DenseVector<ValueType> magneticPermeabilityEMtemp = magneticPermeabilityEM;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp = dielectricPermittivityEM;
    
    magneticPermeabilityEMtemp /= MagneticPermeabilityVacuum;  // calculate the relative magneticPermeabilityEM
    dielectricPermittivityEMtemp /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivityEM
    
    IO::writeVector(magneticPermeabilityEMtemp, filename + ".muEMr", fileFormat);
    IO::writeVector(conductivityEM, filename + ".sigmaEM", fileFormat);
    IO::writeVector(dielectricPermittivityEMtemp, filename + ".epsilonEMr", fileFormat);
    IO::writeVector(porosity, filename + ".porosity", fileFormat);
    IO::writeVector(saturation, filename + ".saturation", fileFormat);    
};

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.EMEM.initializeMatrices")
        HOST_PRINT(comm, "", "Preparation of the Average Matrix\n");

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
void KITGPI::Modelparameter::EMEM<ValueType>::purgeMatrices()
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
void KITGPI::Modelparameter::EMEM<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.EMEM.calculateAveraging")
        this->calculateInverseAveragedEMparameter(magneticPermeabilityEM, inverseMagneticPermeabilityEMAverageYZ, averageMatrixYZ);
        this->calculateInverseAveragedEMparameter(magneticPermeabilityEM, inverseMagneticPermeabilityEMAverageXZ, averageMatrixXZ);
        this->calculateInverseAveragedEMparameter(magneticPermeabilityEM, inverseMagneticPermeabilityEMAverageXY, averageMatrixXY);
        
        this->calculateAveragedEMparameter(conductivityEM, conductivityEMAverageX, averageMatrixX);
        this->calculateAveragedEMparameter(conductivityEM, conductivityEMAverageY, averageMatrixY);
        this->calculateAveragedEMparameter(conductivityEM, conductivityEMAverageZ, averageMatrixZ);
        
        this->calculateAveragedEMparameter(dielectricPermittivityEM, dielectricPermittivityEMAverageX, averageMatrixX);
        this->calculateAveragedEMparameter(dielectricPermittivityEM, dielectricPermittivityEMAverageY, averageMatrixY);
        this->calculateAveragedEMparameter(dielectricPermittivityEM, dielectricPermittivityEMAverageZ, averageMatrixZ);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Get equationType (emem)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::EMEM<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to inverse magneticPermeabilityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getConductivityEMoptical()
{
    COMMON_THROWEXCEPTION("conductivityEMoptical is not set for emem modelling")
    return (conductivityEMoptical);
}

/*! \brief Get reference to dielectricPermittivityEMoptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getDielectricPermittivityEMoptical()
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an emem modelling")
    return (dielectricPermittivityEMoptical);
}

/*! \brief Get reference to tauDisplacementEM
 *
 */
template <typename ValueType> 
ValueType const &KITGPI::Modelparameter::EMEM<ValueType>::getTauDisplacementEM()
{
    COMMON_THROWEXCEPTION("There is no tauDisplacementEM parameter in an emem modelling")
    return (tauDisplacementEM);
}

/*! \brief Get reference to tauDisplacementEM
 *
 */
template <typename ValueType> 
ValueType const &KITGPI::Modelparameter::EMEM<ValueType>::getTauDisplacementEM() const
{
    COMMON_THROWEXCEPTION("There is no tauDisplacementEM parameter in an emem modelling")
    return (tauDisplacementEM);
}

/*! \brief Get reference to tauConductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getTauConductivityEM() const
{
    COMMON_THROWEXCEPTION("There is no tauConductivityEM parameter in an emem modelling")
    return (tauConductivityEM);
}

/*! \brief Get reference to tauDielectricPermittivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getTauDielectricPermittivityEM() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityEM parameter in an emem modelling")
    return (tauDielectricPermittivityEM);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::EMEM<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an emem modelling")
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::EMEM<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an emem modelling")
    return (numRelaxationMechanisms);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> KITGPI::Modelparameter::EMEM<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::EMEM<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::EMEM<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> &KITGPI::Modelparameter::EMEM<ValueType>::operator*=(ValueType const &rhs)
{
    magneticPermeabilityEM *= rhs;
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    porosity *= rhs;
    saturation *= rhs;

    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> KITGPI::Modelparameter::EMEM<ValueType>::operator+(KITGPI::Modelparameter::EMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::EMEM<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> &KITGPI::Modelparameter::EMEM<ValueType>::operator+=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs)
{
    magneticPermeabilityEM += rhs.magneticPermeabilityEM;
    conductivityEM += rhs.conductivityEM;
    dielectricPermittivityEM += rhs.dielectricPermittivityEM;
    porosity += rhs.porosity;
    saturation += rhs.saturation;

    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> KITGPI::Modelparameter::EMEM<ValueType>::operator-(KITGPI::Modelparameter::EMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::EMEM<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> &KITGPI::Modelparameter::EMEM<ValueType>::operator-=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs)
{
    magneticPermeabilityEM -= rhs.magneticPermeabilityEM;
    conductivityEM -= rhs.conductivityEM;
    dielectricPermittivityEM -= rhs.dielectricPermittivityEM;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType> &KITGPI::Modelparameter::EMEM<ValueType>::operator=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs)
{
    magneticPermeabilityEM = rhs.magneticPermeabilityEM;
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    conductivityEMWater = rhs.conductivityEMWater;
    relativeDieletricPeimittivityRockMatrix = rhs.relativeDieletricPeimittivityRockMatrix;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    magneticPermeabilityEM = rhs.getMagneticPermeabilityEM();
    conductivityEM = rhs.getConductivityEM();
    dielectricPermittivityEM = rhs.getDielectricPermittivityEM();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    
    conductivityEMWater = rhs.getConductivityEMWater();
    relativeDieletricPeimittivityRockMatrix = rhs.getRelativeDieletricPeimittivityRockMatrix();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    magneticPermeabilityEM -= rhs.getMagneticPermeabilityEM();
    conductivityEM -= rhs.getConductivityEM();
    dielectricPermittivityEM -= rhs.getDielectricPermittivityEM();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    magneticPermeabilityEM += rhs.getMagneticPermeabilityEM();
    conductivityEM += rhs.getConductivityEM();
    dielectricPermittivityEM += rhs.getDielectricPermittivityEM();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

template class KITGPI::Modelparameter::EMEM<float>;
template class KITGPI::Modelparameter::EMEM<double>;
