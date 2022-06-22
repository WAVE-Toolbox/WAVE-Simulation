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
    /* 9 Parameter in emem modeling: inverseMagneticPermeabilityXY, inverseMagneticPermeabilityYZ, inverseMagneticPermeabilityXY, CaAverageX, CaAverageY, CaAverageZ, CbAverageX, CbAverageY, CbAverageZ */
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
    dielectricPermittivity /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivity
    
    lama::DenseVector<ValueType> mask; //mask to restore vacuum
    mask = dielectricPermittivity-1;
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);
    
    Common::searchAndReplace<ValueType>(electricConductivity, config.get<ValueType>("lowerSigmaTh"), config.get<ValueType>("lowerSigmaTh"), 1);
    Common::searchAndReplace<ValueType>(electricConductivity, config.get<ValueType>("upperSigmaTh"), config.get<ValueType>("upperSigmaTh"), 2);    
    Common::searchAndReplace<ValueType>(dielectricPermittivity, config.get<ValueType>("lowerEpsilonrTh"), config.get<ValueType>("lowerEpsilonrTh"), 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivity, config.get<ValueType>("upperEpsilonrTh"), config.get<ValueType>("upperEpsilonrTh"), 2);
       
    dirtyFlagVelocivityEM = true ;   // If EM-parameters will be changed, velocityEM needs to be redone
    dirtyFlagAveraging = true;      // If EM-parameters will be changed, averaging needs to be redone
      
    if (config.getAndCatch("inversionType", 1) == 3 || config.getAndCatch("parameterisation", 0) == 1 || config.getAndCatch("parameterisation", 0) == 2) {
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("lowerPorosityTh", 0.0), config.getAndCatch("lowerPorosityTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("upperPorosityTh", 1.0), config.getAndCatch("upperPorosityTh", 1.0), 2);

        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("lowerSaturationTh", 0.0), config.getAndCatch("lowerSaturationTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("upperSaturationTh", 1.0), config.getAndCatch("upperSaturationTh", 1.0), 2);
    }
    if (config.getAndCatch("gradientKernel", 0) > 1 && config.getAndCatch("decomposition", 0) == 0) {
        Common::searchAndReplace<ValueType>(reflectivity, config.getAndCatch("lowerReflectivityTh", -1.0), config.getAndCatch("lowerReflectivityTh", -1.0), 1);
        Common::searchAndReplace<ValueType>(reflectivity, config.getAndCatch("upperReflectivityTh", 1.0), config.getAndCatch("upperReflectivityTh", 1.0), 2);
    }
    electricConductivity *= mask;
    porosity *= mask;
    saturation *= mask;
    reflectivity *= mask;
    
    dielectricPermittivity -= 1;
    dielectricPermittivity *= mask;
    dielectricPermittivity += 1;
    dielectricPermittivity *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity   
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);

    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * magneticPermeability;
    modelPerShot.setMagneticPermeability(temp);
    
    temp = shrinkMatrix * dielectricPermittivity;
    modelPerShot.setDielectricPermittivity(temp);
        
    temp = shrinkMatrix * electricConductivity;
    modelPerShot.setElectricConductivity(temp);
    
    temp = shrinkMatrix * porosity;
    modelPerShot.setPorosity(temp);
    
    temp = shrinkMatrix * saturation;
    modelPerShot.setSaturation(temp); 
    
    temp = shrinkMatrix * reflectivity;
    modelPerShot.setReflectivity(temp); 
    
    if (this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        temp = shrinkMatrix * relativeDieletricPeimittivityRockMatrix;
        modelPerShot.setRelativeDieletricPeimittivityRockMatrix(temp);
        
        temp = shrinkMatrix * electricConductivityWater;
        modelPerShot.setElectricConductivityWater(temp);
    }
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
        init(ctx, dist, config.get<ValueType>("mur"), config.get<ValueType>("sigma"), config.get<ValueType>("epsilonr"));
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
    meshingMatrix = lama::zero<lama::CSRSparseMatrix<ValueType>>(variableDist, electricConductivity.getDistributionPtr());
    meshingMatrix.fillFromAssembly(assembly);

    magneticPermeability = meshingMatrix * magneticPermeability;
    electricConductivity = meshingMatrix * electricConductivity;
    dielectricPermittivity = meshingMatrix * dielectricPermittivity;
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeability_const magneticPermeability given as Scalar
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const)
{
    equationType = "emem";
    init(ctx, dist, magneticPermeability_const, electricConductivity_const, dielectricPermittivity_const);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeability_const magneticPermeability given as Scalar
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const)
{
    magneticPermeability_const *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeability
    dielectricPermittivity_const *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity
    
    this->initModelparameter(magneticPermeability, ctx, dist, magneticPermeability_const);
    this->initModelparameter(electricConductivity, ctx, dist, electricConductivity_const);
    this->initModelparameter(dielectricPermittivity, ctx, dist, dielectricPermittivity_const);
    this->initModelparameter(porosity, ctx, dist, 0.0);
    this->initModelparameter(saturation, ctx, dist, 0.0);
    this->initModelparameter(reflectivity, ctx, dist, 0.0);
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
    this->initModelparameter(magneticPermeability, ctx, dist, filename + ".mur", fileFormat);
    this->initModelparameter(electricConductivity, ctx, dist, filename + ".sigma", fileFormat);
    this->initModelparameter(dielectricPermittivity, ctx, dist, filename + ".epsilonr", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        this->initModelparameter(porosity, ctx, dist, filename + ".porosity", fileFormat);
        this->initModelparameter(saturation, ctx, dist, filename + ".saturation", fileFormat);
    } else {
        this->initModelparameter(porosity, ctx, dist, 0.0);
        this->initModelparameter(saturation, ctx, dist, 0.0);
    }
    if (this->getGradientKernel() != 0 && this->getDecomposition() == 0) {
        this->initModelparameter(reflectivity, ctx, dist, filename + ".reflectivity", fileFormat);
    } else {
        this->initModelparameter(reflectivity, ctx, dist, 0.0);
    }
    
    magneticPermeability *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeability
    dielectricPermittivity *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::EMEM<ValueType>::EMEM(const EMEM &rhs)
{
    equationType = rhs.equationType;
    velocivityEM = rhs.velocivityEM;
    magneticPermeability = rhs.magneticPermeability;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    dirtyFlagVelocivityEM = rhs.dirtyFlagVelocivityEM;
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
}

/*! \brief Write model to an external file
 *
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    scai::lama::DenseVector<ValueType> magneticPermeabilitytemp = magneticPermeability;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp = dielectricPermittivity;
    
    magneticPermeabilitytemp /= MagneticPermeabilityVacuum;  // calculate the relative magneticPermeability
    dielectricPermittivitytemp /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivity
    
    IO::writeVector(magneticPermeabilitytemp, filename + ".mur", fileFormat);
    IO::writeVector(electricConductivity, filename + ".sigma", fileFormat);
    IO::writeVector(dielectricPermittivitytemp, filename + ".epsilonr", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        IO::writeVector(porosity, filename + ".porosity", fileFormat);
        IO::writeVector(saturation, filename + ".saturation", fileFormat);
    }
    if (this->getGradientKernel() != 0 && this->getDecomposition() == 0) {
        IO::writeVector(reflectivity, filename + ".reflectivity", fileFormat);
    }
};

//! \brief Initialization of the Averaging matrices
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
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageYZ, averageMatrixYZ);
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageXZ, averageMatrixXZ);
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageXY, averageMatrixXY);
        
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageX, averageMatrixX);
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageY, averageMatrixY);
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageZ, averageMatrixZ);
        
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageX, averageMatrixX);
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageY, averageMatrixY);
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageZ, averageMatrixZ);
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

/*! \brief Get reference to tauElectricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getTauElectricConductivity() const
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivity parameter in an emem modelling")
    return (tauElectricConductivity);
}

/*! \brief Get reference to tauDielectricPermittivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::EMEM<ValueType>::getTauDielectricPermittivity() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivity parameter in an emem modelling")
    return (tauDielectricPermittivity);
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
std::vector<ValueType> KITGPI::Modelparameter::EMEM<ValueType>::getRelaxationFrequency() const
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
    magneticPermeability *= rhs;
    electricConductivity *= rhs;
    dielectricPermittivity *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    reflectivity *= rhs;

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
    magneticPermeability += rhs.magneticPermeability;
    electricConductivity += rhs.electricConductivity;
    dielectricPermittivity += rhs.dielectricPermittivity;
    porosity += rhs.porosity;
    saturation += rhs.saturation;
    reflectivity += rhs.reflectivity;

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
    magneticPermeability -= rhs.magneticPermeability;
    electricConductivity -= rhs.electricConductivity;
    dielectricPermittivity -= rhs.dielectricPermittivity;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;
    reflectivity -= rhs.reflectivity;

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
    magneticPermeability = rhs.magneticPermeability;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    reflectivity = rhs.reflectivity;
    
    electricConductivityWater = rhs.electricConductivityWater;
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
void KITGPI::Modelparameter::EMEM<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability = rhs.getMagneticPermeability();
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    reflectivity = rhs.getReflectivity();
    
    electricConductivityWater = rhs.getElectricConductivityWater();
    relativeDieletricPeimittivityRockMatrix = rhs.getRelativeDieletricPeimittivityRockMatrix();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability -= rhs.getMagneticPermeability();
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    reflectivity -= rhs.getReflectivity();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::EMEM<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability += rhs.getMagneticPermeability();
    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    reflectivity += rhs.getReflectivity();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
}

template class KITGPI::Modelparameter::EMEM<float>;
template class KITGPI::Modelparameter::EMEM<double>;
