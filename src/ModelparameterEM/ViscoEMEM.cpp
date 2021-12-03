#include "ViscoEMEM.hpp"
#include "../IO/IO.hpp"
using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 * 
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::ViscoEMEM<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 13 Parameter in emem modeling: inverseMagneticPermeabilityXY, inverseMagneticPermeabilityYZ, inverseMagneticPermeabilityXZ, dielectricPermittivityAverageX, dielectricPermittivityAverageY, dielectricPermittivityAverageZ, electricConductivityAverageX, electricConductivityAverageY, electricConductivityAverageZ, tauDielectricPermittivityAverageX, tauDielectricPermittivityAverageY, tauDielectricPermittivityAverageZ, tauElectricConductivityAverageX, tauElectricConductivityAverageY, tauElectricConductivityAverageZ */
    IndexType numParameter = 15;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modelparameter for visco-emem modelling
 *
 * Applies Equation 12 from Bohlen 2002 and refreshes the modulus
 *
 \param modelCoordinates coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    //refreshModulus
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
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective = this->getDielectricPermittivityRealEffective();
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective = this->getElectricConductivityRealEffective();
    dielectricPermittivityRealEffective /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivity
    ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * centerFrequencyCPML);
    tauElectricConductivity /= relaxationTime_ref;  // calculate the relative tauElectricConductivity
    
    lama::DenseVector<ValueType> mask; //mask to restore vacuum
    mask = dielectricPermittivityRealEffective-1;
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);
    
    Common::searchAndReplace<ValueType>(electricConductivityRealEffective, config.get<ValueType>("lowerSigmaEMTh"), config.get<ValueType>("lowerSigmaEMTh"), 1);
    Common::searchAndReplace<ValueType>(electricConductivityRealEffective, config.get<ValueType>("upperSigmaEMTh"), config.get<ValueType>("upperSigmaEMTh"), 2);
    Common::searchAndReplace<ValueType>(dielectricPermittivityRealEffective, config.get<ValueType>("lowerEpsilonEMrTh"), config.get<ValueType>("lowerEpsilonEMrTh"), 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityRealEffective, config.get<ValueType>("upperEpsilonEMrTh"), config.get<ValueType>("upperEpsilonEMrTh"), 2);
    Common::searchAndReplace<ValueType>(tauElectricConductivity, config.get<ValueType>("lowerTauSigmaEMrTh"), config.get<ValueType>("lowerTauSigmaEMrTh"), 1);
    Common::searchAndReplace<ValueType>(tauElectricConductivity, config.get<ValueType>("upperTauSigmaEMrTh"), config.get<ValueType>("upperTauSigmaEMrTh"), 2);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivity, config.get<ValueType>("lowerTauEpsilonEMTh"), config.get<ValueType>("lowerTauEpsilonEMTh"), 1);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivity, config.get<ValueType>("upperTauEpsilonEMTh"), config.get<ValueType>("upperTauEpsilonEMTh"), 2);
    
    dirtyFlagAveraging = true;      // If EM-parameters will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;   // the velocity vector is now dirty
    
    if (config.getAndCatch("inversionType", 1) == 3 || config.getAndCatch("parameterisation", 0) == 1 || config.getAndCatch("parameterisation", 0) == 2) {
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("lowerPorosityTh", 0.0), config.getAndCatch("lowerPorosityTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("upperPorosityTh", 1.0), config.getAndCatch("upperPorosityTh", 1.0), 2);

        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("lowerSaturationTh", 0.0), config.getAndCatch("lowerSaturationTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("upperSaturationTh", 1.0), config.getAndCatch("upperSaturationTh", 1.0), 2);
    }
    electricConductivityRealEffective *= mask;
    tauElectricConductivity *= mask;
    tauDielectricPermittivity *= mask;
    porosity *= mask;
    saturation *= mask;
    
    tauElectricConductivity *= relaxationTime_ref;  // calculate the real tauElectricConductivity
    dielectricPermittivityRealEffective -= 1;
    dielectricPermittivityRealEffective *= mask;
    dielectricPermittivityRealEffective += 1;
    dielectricPermittivityRealEffective *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity
    
    IndexType calculateType = 2;
    dielectricPermittivity = this->getDielectricPermittivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType);
    electricConductivity = this->getElectricConductivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType);
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);

    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * dielectricPermittivity;
    modelPerShot.setDielectricPermittivity(temp);
        
    temp = shrinkMatrix * electricConductivity;
    modelPerShot.setElectricConductivity(temp);
    
    temp = shrinkMatrix * tauDielectricPermittivity;
    modelPerShot.setTauDielectricPermittivity(temp);
        
    temp = shrinkMatrix * tauElectricConductivity;
    modelPerShot.setTauElectricConductivity(temp);
    
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
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivity.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivity().getDistributionPtr();
//     auto comm = dist.getCommunicatorPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate, boundaryWidth);
    scai::lama::SparseVector<ValueType> restoreVector;
    restoreVector = 1.0 - eraseVector;
    
    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * modelPerShot.getDielectricPermittivity(); //transform pershot into big model
    temp *= restoreVector;
    dielectricPermittivity *= eraseVector;
    dielectricPermittivity += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getElectricConductivity(); //transform pershot into big model
    temp *= restoreVector;
    electricConductivity *= eraseVector;
    electricConductivity += temp; //take over the values
    
    temp = shrinkMatrix * modelPerShot.getTauDielectricPermittivity(); //transform pershot into big model
    temp *= restoreVector;
    tauDielectricPermittivity *= eraseVector;
    tauDielectricPermittivity += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getTauElectricConductivity(); //transform pershot into big model
    temp *= restoreVector;
    tauElectricConductivity *= eraseVector;
    tauElectricConductivity += temp; //take over the values
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType>::ViscoEMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "viscoemem";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2, "Read variable model not available for ViscoEMEM, variable grid is not available here!")
    
    IndexType numRelaxationMechanisms_in = config.get<IndexType>("numRelaxationMechanisms");
    SCAI_ASSERT(numRelaxationMechanisms_in <= 4, "numRelaxationMechanisms more than 4 is not available here!")
    std::vector<ValueType> relaxationFrequency_in(numRelaxationMechanisms_in, 0);
    for (int l=0; l<numRelaxationMechanisms_in; l++) {
        if (l==0)
            relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency");
        if (l==1)
            relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency2");
        if (l==2)
            relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency3");
        if (l==3)
            relaxationFrequency_in[l] = config.get<ValueType>("relaxationFrequency4");
    }
    if (config.get<IndexType>("ModelRead") == 1) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (viscoemem) parameter from file...\n");

        initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML")); 
        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("muEMr"), config.get<ValueType>("sigmaEM"), config.get<ValueType>("epsilonEMr"), config.get<ValueType>("tauSigmaEMr"), config.get<ValueType>("tauEpsilonEM"), numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeability_const magneticPermeability given as Scalar
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 \param tauElectricConductivity_const tauElectricConductivity given as Scalar
 \param tauDielectricPermittivity_const TtauDielectricPermittivity given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType>::ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    equationType = "viscoemem";
    init(ctx, dist, magneticPermeability_const, electricConductivity_const, dielectricPermittivity_const, tauElectricConductivity_const, tauDielectricPermittivity_const, numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the five given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeability_const magneticPermeability given as Scalar
 \param electricConductivity_const electricConductivity given as Scalar
 \param dielectricPermittivity_const dielectricPermittivity given as Scalar
 \param tauElectricConductivity_const tauElectricConductivity given as Scalar
 \param tauDielectricPermittivity_const TtauDielectricPermittivity given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
    
    dielectricPermittivity_const *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity
    magneticPermeability_const *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeability
    
    this->initModelparameter(magneticPermeability, ctx, dist, magneticPermeability_const);
    this->initModelparameter(electricConductivity, ctx, dist, electricConductivity_const);
    this->initModelparameter(dielectricPermittivity, ctx, dist, dielectricPermittivity_const);
    
    ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * centerFrequencyCPML);
    tauElectricConductivity_const *= relaxationTime_ref; // calculate the real tauElectricConductivity
    
    this->initModelparameter(tauElectricConductivity, ctx, dist, tauElectricConductivity_const);
    this->initModelparameter(tauDielectricPermittivity, ctx, dist, tauDielectricPermittivity_const);
    this->initModelparameter(porosity, ctx, dist, 0.0);
    this->initModelparameter(saturation, ctx, dist, 0.0);
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
KITGPI::Modelparameter::ViscoEMEM<ValueType>::ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "viscoemem";
    init(ctx, dist, filename, fileFormat);
}

/*! \brief Initialisator that is reading Model-Vector from an external files and calculates electricConductivity and dielectricPermittivity
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective;
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective;
    
    this->initModelparameter(magneticPermeability, ctx, dist, filename + ".muEMr", fileFormat);
    this->initModelparameter(electricConductivityRealEffective, ctx, dist, filename + ".sigmaEM", fileFormat);
    this->initModelparameter(dielectricPermittivityRealEffective, ctx, dist, filename + ".epsilonEMr", fileFormat);
    this->initModelparameter(tauElectricConductivity, ctx, dist, filename + ".tauSigmaEMr", fileFormat);
    this->initModelparameter(tauDielectricPermittivity, ctx, dist, filename + ".tauEpsilonEM", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        this->initModelparameter(porosity, ctx, dist, filename + ".porosity", fileFormat);
        this->initModelparameter(saturation, ctx, dist, filename + ".saturation", fileFormat);
    } else {
        this->initModelparameter(porosity, ctx, dist, 0.0);
        this->initModelparameter(saturation, ctx, dist, 0.0);
    }
            
    magneticPermeability *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeability
    ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * centerFrequencyCPML);
    tauElectricConductivity *= relaxationTime_ref; // calculate the real tauElectricConductivity
    dielectricPermittivityRealEffective *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivity
    
    IndexType calculateType = 2;
    dielectricPermittivity = this->getDielectricPermittivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType);
    electricConductivity = this->getElectricConductivityStatic(dielectricPermittivityRealEffective, electricConductivityRealEffective, calculateType);
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType>::ViscoEMEM(const ViscoEMEM &rhs)
{
    equationType = rhs.equationType;
    velocivityEM = rhs.velocivityEM;
    magneticPermeability = rhs.magneticPermeability;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;   
    dirtyFlagVelocivityEM = rhs.dirtyFlagVelocivityEM;
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    tauElectricConductivity = rhs.tauElectricConductivity;
    tauDielectricPermittivity = rhs.tauDielectricPermittivity;
    
    relaxationFrequency = rhs.relaxationFrequency;    
    numRelaxationMechanisms = rhs.numRelaxationMechanisms; 
    centerFrequencyCPML = rhs.centerFrequencyCPML;
}

/*! \brief Write model to an external file
 *base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    scai::lama::DenseVector<ValueType> magneticPermeabilitytemp = magneticPermeability;
    scai::lama::DenseVector<ValueType> tauElectricConductivitytemp = tauElectricConductivity;
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective = this->getDielectricPermittivityRealEffective();
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective = this->getElectricConductivityRealEffective();
    
    magneticPermeabilitytemp /= MagneticPermeabilityVacuum;  // calculate the relative magneticPermeability
    ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * centerFrequencyCPML);
    tauElectricConductivitytemp /= relaxationTime_ref;  // calculate the relative tauElectricConductivity
    dielectricPermittivityRealEffective /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivity
    
    IO::writeVector(magneticPermeabilitytemp, filename + ".muEMr", fileFormat);
    IO::writeVector(electricConductivityRealEffective, filename + ".sigmaEM", fileFormat);
    IO::writeVector(dielectricPermittivityRealEffective, filename + ".epsilonEMr", fileFormat);
    IO::writeVector(tauElectricConductivitytemp, filename + ".tauSigmaEMr", fileFormat);
    IO::writeVector(tauDielectricPermittivity, filename + ".tauEpsilonEM", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        IO::writeVector(porosity, filename + ".porosity", fileFormat);
        IO::writeVector(saturation, filename + ".saturation", fileFormat);
    }
};

//! \brief Initialization of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates coordinate class object
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoEMEM.initializeMatrices")
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
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::purgeMatrices()
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
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoEMEM.calculateAveraging")
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageYZ, averageMatrixYZ);
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageXZ, averageMatrixXZ);
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageXY, averageMatrixXY);
                
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageX, averageMatrixX);
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageY, averageMatrixY);
        this->calcAveragedParameter(electricConductivity, electricConductivityAverageZ, averageMatrixZ);
        
        this->calcAveragedParameter(tauElectricConductivity, tauElectricConductivityAverageX, averageMatrixX);
        this->calcAveragedParameter(tauElectricConductivity, tauElectricConductivityAverageY, averageMatrixY);
        this->calcAveragedParameter(tauElectricConductivity, tauElectricConductivityAverageZ, averageMatrixZ);
        
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageX, averageMatrixX);
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageY, averageMatrixY);
        this->calcAveragedParameter(dielectricPermittivity, dielectricPermittivityAverageZ, averageMatrixZ);
        
        this->calcAveragedParameter(tauDielectricPermittivity, tauDielectricPermittivityAverageX, averageMatrixX);
        this->calcAveragedParameter(tauDielectricPermittivity, tauDielectricPermittivityAverageY, averageMatrixY);
        this->calcAveragedParameter(tauDielectricPermittivity, tauDielectricPermittivityAverageZ, averageMatrixZ);
        
        dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-emem simulation")
    }
    if (relaxationFrequency_in[0] <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-emem simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
    centerFrequencyCPML = centerFrequencyCPML_in;
}

/*! \brief Get equationType (viscoemem)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::ViscoEMEM<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get const reference to EM-wave velocity
 * 
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoEMEM<ValueType>::getVelocityEM()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagVelocivityEM) {
        HOST_PRINT(dielectricPermittivity.getDistributionPtr()->getCommunicatorPtr(), "", "EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivityRealEffective\n");
        scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective = this->getDielectricPermittivityRealEffective();
        this->calcVelocityFromModulus(dielectricPermittivityRealEffective, magneticPermeability, velocivityEM);
        dirtyFlagVelocivityEM = false;
    }
    return (velocivityEM);
}

/*! \brief Get const reference to EM-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoEMEM<ValueType>::getVelocityEM() const
{
    SCAI_ASSERT(dirtyFlagVelocivityEM == false, "EM-Wave velocity has to be recalculated! ");
    return (velocivityEM);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::ViscoEMEM<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::ViscoEMEM<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> &KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator*=(ValueType const &rhs)
{
    magneticPermeability *= rhs;
    electricConductivity *= rhs;
    dielectricPermittivity *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity *= rhs;
    tauDielectricPermittivity *= rhs;
    
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator+(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoEMEM<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> &KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator+=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs)
{
    magneticPermeability += rhs.magneticPermeability;
    electricConductivity += rhs.electricConductivity;
    dielectricPermittivity += rhs.dielectricPermittivity;
    porosity += rhs.porosity;
    saturation += rhs.saturation;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity += rhs.tauElectricConductivity;
    tauDielectricPermittivity += rhs.tauDielectricPermittivity;
    
    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator-(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoEMEM<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> &KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator-=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs)
{
    magneticPermeability -= rhs.magneticPermeability;
    electricConductivity -= rhs.electricConductivity;
    dielectricPermittivity -= rhs.dielectricPermittivity;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity -= rhs.tauElectricConductivity;
    tauDielectricPermittivity -= rhs.tauDielectricPermittivity;
    
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoEMEM<ValueType> &KITGPI::Modelparameter::ViscoEMEM<ValueType>::operator=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs)
{
    magneticPermeability = rhs.magneticPermeability;
    electricConductivity = rhs.electricConductivity;
    dielectricPermittivity = rhs.dielectricPermittivity;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    electricConductivityWater = rhs.electricConductivityWater;
    relativeDieletricPeimittivityRockMatrix = rhs.relativeDieletricPeimittivityRockMatrix;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity = rhs.tauElectricConductivity;
    tauDielectricPermittivity = rhs.tauDielectricPermittivity;
    
    relaxationFrequency = rhs.relaxationFrequency;    
    numRelaxationMechanisms = rhs.numRelaxationMechanisms; 
    centerFrequencyCPML = rhs.centerFrequencyCPML;
    
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability = rhs.getMagneticPermeability();
    electricConductivity = rhs.getElectricConductivity();
    dielectricPermittivity = rhs.getDielectricPermittivity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    
    electricConductivityWater = rhs.getElectricConductivityWater();
    relativeDieletricPeimittivityRockMatrix = rhs.getRelativeDieletricPeimittivityRockMatrix();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity = rhs.getTauElectricConductivity();
    tauDielectricPermittivity = rhs.getTauDielectricPermittivity();
    relaxationFrequency = rhs.getRelaxationFrequency();
    numRelaxationMechanisms = rhs.getNumRelaxationMechanisms(); 
    centerFrequencyCPML = rhs.getCenterFrequencyCPML();
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability -= rhs.getMagneticPermeability();
    electricConductivity -= rhs.getElectricConductivity();
    dielectricPermittivity -= rhs.getDielectricPermittivity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity -= rhs.getTauElectricConductivity();
    tauDielectricPermittivity -= rhs.getTauDielectricPermittivity();
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoEMEM<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    magneticPermeability += rhs.getMagneticPermeability();
    electricConductivity += rhs.getElectricConductivity();
    dielectricPermittivity += rhs.getDielectricPermittivity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauElectricConductivity += rhs.getTauElectricConductivity();
    tauDielectricPermittivity += rhs.getTauDielectricPermittivity();
}

template class KITGPI::Modelparameter::ViscoEMEM<float>;
template class KITGPI::Modelparameter::ViscoEMEM<double>;
