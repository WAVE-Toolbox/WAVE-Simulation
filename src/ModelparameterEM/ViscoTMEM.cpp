#include "ViscoTMEM.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 * 
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::ViscoTMEM<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 6 Parameter in tmem modeling: inverseMagneticPermeabilityYZ, inverseMagneticPermeabilityXY, dielectricPermittivity, electricConductivity, tauDielectricPermittivity, tauElectricConductivity */
    IndexType numParameter = 6;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modelparameter for modelling
 *
 * Refreshes the modulus and calculates average values on staggered grid
 *
 \param modelCoordinates coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::applyThresholds(Configuration::Configuration const &config)
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
    
    Common::searchAndReplace<ValueType>(electricConductivityRealEffective, config.get<ValueType>("lowerSigmaTh"), config.get<ValueType>("lowerSigmaTh"), 1);
    Common::searchAndReplace<ValueType>(electricConductivityRealEffective, config.get<ValueType>("upperSigmaTh"), config.get<ValueType>("upperSigmaTh"), 2);
    Common::searchAndReplace<ValueType>(dielectricPermittivityRealEffective, config.get<ValueType>("lowerEpsilonrTh"), config.get<ValueType>("lowerEpsilonrTh"), 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityRealEffective, config.get<ValueType>("upperEpsilonrTh"), config.get<ValueType>("upperEpsilonrTh"), 2);
    Common::searchAndReplace<ValueType>(tauElectricConductivity, config.get<ValueType>("lowerTauSigmarTh"), config.get<ValueType>("lowerTauSigmarTh"), 1);
    Common::searchAndReplace<ValueType>(tauElectricConductivity, config.get<ValueType>("upperTauSigmarTh"), config.get<ValueType>("upperTauSigmarTh"), 2);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivity, config.get<ValueType>("lowerTauEpsilonTh"), config.get<ValueType>("lowerTauEpsilonTh"), 1);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivity, config.get<ValueType>("upperTauEpsilonTh"), config.get<ValueType>("upperTauEpsilonTh"), 2);
    
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
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
    
    temp = shrinkMatrix * tauDielectricPermittivity;
    modelPerShot.setTauDielectricPermittivity(temp);
        
    temp = shrinkMatrix * tauElectricConductivity;
    modelPerShot.setTauElectricConductivity(temp);
    
    temp = shrinkMatrix * porosity;
    modelPerShot.setPorosity(temp);
    
    temp = shrinkMatrix * saturation;
    modelPerShot.setSaturation(temp);
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "viscotmem";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{    
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2 || config.get<IndexType>("UseVariableGrid") == 1, "Read variable model (ModelRead=2) not available if regular grid is chosen!")

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
    if ((config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 0) || (config.get<IndexType>("ModelRead") == 2)) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (viscotmem) parameter from file...\n");

        initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML")); 
        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else if (config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 1) {

        Acquisition::Coordinates<ValueType> regularCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        dmemo::DistributionPtr regularDist(new dmemo::BlockDistribution(regularCoordinates.getNGridpoints(), dist->getCommunicatorPtr()));

        initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML")); 
        init(ctx, regularDist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "reading regular model finished\n\n")

        init(dist, modelCoordinates, regularCoordinates);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "initialising model on discontineous grid finished\n")

    } else {
        init(ctx, dist, config.get<ValueType>("mur"), config.get<ValueType>("sigma"), config.get<ValueType>("epsilonr"), config.get<ValueType>("tauSigmar"), config.get<ValueType>("tauEpsilon"), numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));
    }
}

/*! \brief initialisation function which creates a variable grid model on top of a regular model
    \param model regular input model
    \param variableDist Distribution for a variable grid
    \param variableCoordinates Coordinate Class of a Variable Grid
    \param regularCoordinates Coordinate Class of a regular Grid
    */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates)
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
 \param tauElectricConductivity_const tauElectricConductivity given as Scalar
 \param tauDielectricPermittivity_const TtauDielectricPermittivity given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    equationType = "viscotmem";
    init(ctx, dist, magneticPermeability_const, electricConductivity_const, dielectricPermittivity_const, tauElectricConductivity_const, tauDielectricPermittivity_const, numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
}

/*! \brief Initialisation that is generating a homogeneous model
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
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
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "viscotmem";
    init(ctx, dist, filename, fileFormat);
}

/*! \brief Initialisator that is reading Model-Vector from an external files and calculates electricConductivity and dielectricPermittivity
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective;
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective;
    
    this->initModelparameter(magneticPermeability, ctx, dist, filename + ".mur", fileFormat);
    this->initModelparameter(electricConductivityRealEffective, ctx, dist, filename + ".sigma", fileFormat);
    this->initModelparameter(dielectricPermittivityRealEffective, ctx, dist, filename + ".epsilonr", fileFormat);
    this->initModelparameter(tauElectricConductivity, ctx, dist, filename + ".tauSigmar", fileFormat);
    this->initModelparameter(tauDielectricPermittivity, ctx, dist, filename + ".tauEpsilon", fileFormat);
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
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(const ViscoTMEM &rhs)
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
 *
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    scai::lama::DenseVector<ValueType> magneticPermeabilitytemp = magneticPermeability;
    scai::lama::DenseVector<ValueType> tauElectricConductivitytemp = tauElectricConductivity;
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective = this->getDielectricPermittivityRealEffective();
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective = this->getElectricConductivityRealEffective();
    
    magneticPermeabilitytemp /= MagneticPermeabilityVacuum;  // calculate the relative magneticPermeability
    ValueType relaxationTime_ref = 1.0 / (2.0 * M_PI * centerFrequencyCPML);
    tauElectricConductivitytemp /= relaxationTime_ref;  // calculate the relative tauElectricConductivity
    dielectricPermittivityRealEffective /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivity
    
    IO::writeVector(magneticPermeabilitytemp, filename + ".mur", fileFormat);
    IO::writeVector(electricConductivityRealEffective, filename + ".sigma", fileFormat);
    IO::writeVector(dielectricPermittivityRealEffective, filename + ".epsilonr", fileFormat);
    IO::writeVector(tauElectricConductivitytemp, filename + ".tauSigmar", fileFormat);
    IO::writeVector(tauDielectricPermittivity, filename + ".tauEpsilon", fileFormat);
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
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoTMEM.initializeMatrices")
        HOST_PRINT(comm, "", "Preparation of the Average Matrix\n");

        this->calcAverageMatrixX(modelCoordinates, dist);
        this->calcAverageMatrixY(modelCoordinates, dist);

        averageMatrixX.setContextPtr(ctx);
        averageMatrixY.setContextPtr(ctx);
    }
}

//! \brief Purge Averaging matrices to free memory
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::purgeMatrices()
{
    averageMatrixX.purge();
    averageMatrixY.purge();
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoTMEM.calculateAveraging")
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageXZ, averageMatrixX);
        this->calcInverseAveragedParameter(magneticPermeability, inverseMagneticPermeabilityAverageYZ, averageMatrixY);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-tmem simulation")
    }
    if (relaxationFrequency_in[0] <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-tmem simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
    centerFrequencyCPML = centerFrequencyCPML_in;
}

/*! \brief Get equationType (viscotmem)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::ViscoTMEM<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get const reference to EM-wave velocity
 * 
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoTMEM<ValueType>::getVelocityEM()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoTMEM<ValueType>::getVelocityEM() const
{
    SCAI_ASSERT(dirtyFlagVelocivityEM == false, "EM-Wave velocity has to be recalculated! ");
    return (velocivityEM);
}

/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::ViscoTMEM<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::ViscoTMEM<ValueType> rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> &KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator*=(ValueType const &rhs)
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
KITGPI::Modelparameter::ViscoTMEM<ValueType> KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator+(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoTMEM<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> &KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator+=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
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
KITGPI::Modelparameter::ViscoTMEM<ValueType> KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator-(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoTMEM<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> &KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator-=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
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
KITGPI::Modelparameter::ViscoTMEM<ValueType> &KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
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

template class KITGPI::Modelparameter::ViscoTMEM<float>;
template class KITGPI::Modelparameter::ViscoTMEM<double>;
