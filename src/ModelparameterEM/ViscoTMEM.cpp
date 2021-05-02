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
    /* 6 Parameter in emem modeling: inverseMagneticPermeabilityEMYZ, inverseMagneticPermeabilityEMXY, CaAverageZ, CbAverageZ, Cc, CdAverageZ */
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
    this->getTauDisplacementEM();
    this->getConductivityEMoptical();
    this->getDielectricPermittivityEMoptical();
    this->getVelocityEM();
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    calculateAveraging();
    purgeMatrices();
    HOST_PRINT(comm, "", "Model ready!\n\n");

//     std::cout << "prepareForModelling relaxationFrequency = " << relaxationFrequency << "\n" << std::endl;
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    dielectricPermittivityEM /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivityEM
    tauConductivityEM /= this->getTauDisplacementEM();  // calculate the relative tauConductivityEM
    
    lama::DenseVector<ValueType> mask; //mask to restore vacuum
    mask = dielectricPermittivityEM-1;
    mask.unaryOp(mask, common::UnaryOp::SIGN);
    mask.unaryOp(mask, common::UnaryOp::ABS);
    
    Common::searchAndReplace<ValueType>(conductivityEM, config.get<ValueType>("lowerSigmaEMTh"), config.get<ValueType>("lowerSigmaEMTh"), 1);
    Common::searchAndReplace<ValueType>(conductivityEM, config.get<ValueType>("upperSigmaEMTh"), config.get<ValueType>("upperSigmaEMTh"), 2);
    Common::searchAndReplace<ValueType>(dielectricPermittivityEM, config.get<ValueType>("lowerEpsilonEMrTh"), config.get<ValueType>("lowerEpsilonEMrTh"), 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityEM, config.get<ValueType>("upperEpsilonEMrTh"), config.get<ValueType>("upperEpsilonEMrTh"), 2);
    Common::searchAndReplace<ValueType>(tauConductivityEM, config.get<ValueType>("lowerTauSigmaEMrTh"), config.get<ValueType>("lowerTauSigmaEMrTh"), 1);
    Common::searchAndReplace<ValueType>(tauConductivityEM, config.get<ValueType>("upperTauSigmaEMrTh"), config.get<ValueType>("upperTauSigmaEMrTh"), 2);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivityEM, config.get<ValueType>("lowerTauEpsilonEMTh"), config.get<ValueType>("lowerTauEpsilonEMTh"), 1);
    Common::searchAndReplace<ValueType>(tauDielectricPermittivityEM, config.get<ValueType>("upperTauEpsilonEMTh"), config.get<ValueType>("upperTauEpsilonEMTh"), 2);
    
    dirtyFlagAveraging = true;      // If EM-parameters will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;   // the velocity vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the conductivityEM vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the dielectricPermittivityEMoptical vector is now dirty
          
    if (config.get<IndexType>("inversionType") == 3 || config.get<IndexType>("parameterisation") == 1 || config.get<IndexType>("parameterisation") == 2) {
        Common::searchAndReplace<ValueType>(porosity, config.get<ValueType>("lowerPorosityTh"), config.get<ValueType>("lowerPorosityTh"), 1);
        Common::searchAndReplace<ValueType>(porosity, config.get<ValueType>("upperPorosityTh"), config.get<ValueType>("upperPorosityTh"), 2);

        Common::searchAndReplace<ValueType>(saturation, config.get<ValueType>("lowerSaturationTh"), config.get<ValueType>("lowerSaturationTh"), 1);
        Common::searchAndReplace<ValueType>(saturation, config.get<ValueType>("upperSaturationTh"), config.get<ValueType>("upperSaturationTh"), 2);
    }
    conductivityEM *= mask;
    tauConductivityEM *= mask;
    tauDielectricPermittivityEM *= mask;
    porosity *= mask;
    saturation *= mask;
    
    dielectricPermittivityEM -= 1;
    dielectricPermittivityEM *= mask;
    dielectricPermittivityEM += 1;
    dielectricPermittivityEM *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM
    tauConductivityEM *= this->getTauDisplacementEM();  // calculate the real tauConductivityEM
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivityEM().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);

    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * dielectricPermittivityEM;
    modelPerShot.setDielectricPermittivityEM(temp);
        
    temp = shrinkMatrix * conductivityEM;
    modelPerShot.setConductivityEM(temp);
    
    temp = shrinkMatrix * tauDielectricPermittivityEM;
    modelPerShot.setTauDielectricPermittivityEM(temp);
        
    temp = shrinkMatrix * tauConductivityEM;
    modelPerShot.setTauConductivityEM(temp);
    
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    auto distBig = dielectricPermittivityEM.getDistributionPtr();
    auto dist = modelPerShot.getDielectricPermittivityEM().getDistributionPtr();
//     auto comm = dist.getCommunicatorPtr();

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
    
    temp = shrinkMatrix * modelPerShot.getTauDielectricPermittivityEM(); //transform pershot into big model
    temp *= restoreVector;
    tauDielectricPermittivityEM *= eraseVector;
    tauDielectricPermittivityEM += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getTauConductivityEM(); //transform pershot into big model
    temp *= restoreVector;
    tauConductivityEM *= eraseVector;
    tauConductivityEM += temp; //take over the values
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

    if ((config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 0) || (config.get<IndexType>("ModelRead") == 2)) {

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (viscotmem) parameter from file...\n");

        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency")); 
        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else if (config.get<IndexType>("ModelRead") == 1 && config.get<IndexType>("UseVariableGrid") == 1) {

        Acquisition::Coordinates<ValueType> regularCoordinates(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        dmemo::DistributionPtr regularDist(new dmemo::BlockDistribution(regularCoordinates.getNGridpoints(), dist->getCommunicatorPtr()));

        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency")); 
        init(ctx, regularDist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "reading regular model finished\n\n")

        init(dist, modelCoordinates, regularCoordinates);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "initialising model on discontineous grid finished\n")

    } else {
        init(ctx, dist, config.get<ValueType>("muEMr"), config.get<ValueType>("sigmaEM"), config.get<ValueType>("epsilonEMr"), config.get<ValueType>("tauSigmaEMr"), config.get<ValueType>("tauEpsilonEM"), config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));
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
 \param tauConductivityEM_const tauConductivityEM given as Scalar
 \param tauDielectricPermittivityEM_const TtauDielectricPermittivityEM given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType tauConductivityEM_const, ValueType tauDielectricPermittivityEM_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    equationType = "viscotmem";
    init(ctx, dist, magneticPermeabilityEM_const, conductivityEM_const, dielectricPermittivityEM_const, tauConductivityEM_const, tauDielectricPermittivityEM_const, numRelaxationMechanisms_in, relaxationFrequency_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param magneticPermeabilityEM_const magneticPermeabilityEM given as Scalar
 \param conductivityEM_const conductivityEM given as Scalar
 \param dielectricPermittivityEM_const dielectricPermittivityEM given as Scalar
 \param tauConductivityEM_const tauConductivityEM given as Scalar
 \param tauDielectricPermittivityEM_const TtauDielectricPermittivityEM given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType tauConductivityEM_const, ValueType tauDielectricPermittivityEM_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    dielectricPermittivityEM_const *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM
    magneticPermeabilityEM_const *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeabilityEM
    
    this->initModelparameter(magneticPermeabilityEM, ctx, dist, magneticPermeabilityEM_const);
    this->initModelparameter(conductivityEM, ctx, dist, conductivityEM_const);
    this->initModelparameter(dielectricPermittivityEM, ctx, dist, dielectricPermittivityEM_const);
    
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    tauConductivityEM_const *= this->getTauDisplacementEM(); // calculate the real tauConductivityEM
    
    this->initModelparameter(tauConductivityEM, ctx, dist, tauConductivityEM_const);
    this->initModelparameter(tauDielectricPermittivityEM, ctx, dist, tauDielectricPermittivityEM_const);
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
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    this->initModelparameter(magneticPermeabilityEM, ctx, dist, filename + ".muEMr", fileFormat);
    this->initModelparameter(conductivityEM, ctx, dist, filename + ".sigmaEM", fileFormat);
    this->initModelparameter(dielectricPermittivityEM, ctx, dist, filename + ".epsilonEMr", fileFormat);
    this->initModelparameter(tauConductivityEM, ctx, dist, filename + ".tauSigmaEMr", fileFormat);
    this->initModelparameter(tauDielectricPermittivityEM, ctx, dist, filename + ".tauEpsilonEM", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        this->initModelparameter(porosity, ctx, dist, filename + ".porosity", fileFormat);
        this->initModelparameter(saturation, ctx, dist, filename + ".saturation", fileFormat);
    } else {
        this->initModelparameter(porosity, ctx, dist, 0.0);
        this->initModelparameter(saturation, ctx, dist, 0.0);
    }
            
    magneticPermeabilityEM *= MagneticPermeabilityVacuum;  // calculate the real magneticPermeabilityEM
    dielectricPermittivityEM *= DielectricPermittivityVacuum;  // calculate the real dielectricPermittivityEM
    tauConductivityEM *= this->getTauDisplacementEM(); // calculate the real tauConductivityEM
}

//! \brief Copy constructor
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType>::ViscoTMEM(const ViscoTMEM &rhs)
{
    equationType = rhs.equationType;
    velocivityEM = rhs.velocivityEM;
    magneticPermeabilityEM = rhs.magneticPermeabilityEM;
    conductivityEM = rhs.conductivityEM;
    dielectricPermittivityEM = rhs.dielectricPermittivityEM;
    dirtyFlagVelocivityEM = rhs.dirtyFlagVelocivityEM;
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    conductivityEMoptical = rhs.conductivityEMoptical;
    dielectricPermittivityEMoptical = rhs.dielectricPermittivityEMoptical;
    tauConductivityEM = rhs.tauConductivityEM;
    tauDielectricPermittivityEM = rhs.tauDielectricPermittivityEM;
    tauDisplacementEM = rhs.tauDisplacementEM;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    dirtyFlagConductivityEMoptical = rhs.dirtyFlagConductivityEMoptical;
    dirtyFlagDielectricPermittivityEMoptical = rhs.dirtyFlagDielectricPermittivityEMoptical;
}

/*! \brief Write model to an external file
 *
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    scai::lama::DenseVector<ValueType> magneticPermeabilityEMtemp = magneticPermeabilityEM;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp = dielectricPermittivityEM;
    scai::lama::DenseVector<ValueType> tauConductivityEMtemp = tauConductivityEM;
    
    magneticPermeabilityEMtemp /= MagneticPermeabilityVacuum;  // calculate the relative magneticPermeabilityEM
    dielectricPermittivityEMtemp /= DielectricPermittivityVacuum;  // calculate the relative dielectricPermittivityEM
    tauConductivityEMtemp /= this->getTauDisplacementEM();  // calculate the relative tauConductivityEM
    
    IO::writeVector(magneticPermeabilityEMtemp, filename + ".muEMr", fileFormat);
    IO::writeVector(conductivityEM, filename + ".sigmaEM", fileFormat);
    IO::writeVector(dielectricPermittivityEMtemp, filename + ".epsilonEMr", fileFormat);
    IO::writeVector(tauConductivityEMtemp, filename + ".tauSigmaEMr", fileFormat);
    IO::writeVector(tauDielectricPermittivityEM, filename + ".tauEpsilonEM", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        IO::writeVector(porosity, filename + ".porosity", fileFormat);
        IO::writeVector(saturation, filename + ".saturation", fileFormat);
    }
};

//! \brief Initializsation of the Averaging matrices
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
        this->calculateInverseAveragedEMparameter(magneticPermeabilityEM, inverseMagneticPermeabilityEMAverageXZ, averageMatrixX);
        this->calculateInverseAveragedEMparameter(magneticPermeabilityEM, inverseMagneticPermeabilityEMAverageYZ, averageMatrixY);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-emem simulation")
    }
    if (relaxationFrequency_in <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-emem simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
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
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeabilityEM and dielectricPermittivityEMoptical.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoTMEM<ValueType>::getVelocityEM()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagVelocivityEM) {
        HOST_PRINT(dielectricPermittivityEMoptical.getDistributionPtr()->getCommunicatorPtr(), "", "EM-Wave velocity will be calculated from magneticPermeabilityEM and dielectricPermittivityEMoptical\n");
        this->calcVelocityFromModulus(dielectricPermittivityEMoptical, magneticPermeabilityEM, velocivityEM);
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
    magneticPermeabilityEM *= rhs;
    conductivityEM *= rhs;
    dielectricPermittivityEM *= rhs;
    porosity *= rhs;
    saturation *= rhs;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauConductivityEM *= rhs;
    tauDielectricPermittivityEM *= rhs;
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
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
    magneticPermeabilityEM += rhs.magneticPermeabilityEM;
    conductivityEM += rhs.conductivityEM;
    dielectricPermittivityEM += rhs.dielectricPermittivityEM;
    porosity += rhs.porosity;
    saturation += rhs.saturation;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauConductivityEM += rhs.tauConductivityEM;
    tauDielectricPermittivityEM += rhs.tauDielectricPermittivityEM;
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
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
    magneticPermeabilityEM -= rhs.magneticPermeabilityEM;
    conductivityEM -= rhs.conductivityEM;
    dielectricPermittivityEM -= rhs.dielectricPermittivityEM;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauConductivityEM -= rhs.tauConductivityEM;
    tauDielectricPermittivityEM -= rhs.tauDielectricPermittivityEM;
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoTMEM<ValueType> &KITGPI::Modelparameter::ViscoTMEM<ValueType>::operator=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs)
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
    
    tauConductivityEM = rhs.tauConductivityEM;
    tauDielectricPermittivityEM = rhs.tauDielectricPermittivityEM;
    tauDisplacementEM = rhs.tauDisplacementEM;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
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
    
    tauConductivityEM = rhs.getTauConductivityEM();
    tauDielectricPermittivityEM = rhs.getTauDielectricPermittivityEM();
    tauDisplacementEM = rhs.getTauDisplacementEM();
    relaxationFrequency = rhs.getRelaxationFrequency();
    numRelaxationMechanisms = rhs.getNumRelaxationMechanisms();
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    magneticPermeabilityEM -= rhs.getMagneticPermeabilityEM();
    conductivityEM -= rhs.getConductivityEM();
    dielectricPermittivityEM -= rhs.getDielectricPermittivityEM();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauConductivityEM -= rhs.getTauConductivityEM();
    tauDielectricPermittivityEM -= rhs.getTauDielectricPermittivityEM();
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoTMEM<ValueType>::plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    magneticPermeabilityEM += rhs.getMagneticPermeabilityEM();
    conductivityEM += rhs.getConductivityEM();
    dielectricPermittivityEM += rhs.getDielectricPermittivityEM();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagAveraging = true;
    dirtyFlagVelocivityEM = true;
    
    tauConductivityEM += rhs.getTauConductivityEM();
    tauDielectricPermittivityEM += rhs.getTauDielectricPermittivityEM();
    
    dirtyFlagConductivityEMoptical = true;
    dirtyFlagDielectricPermittivityEMoptical = true;
}

template class KITGPI::Modelparameter::ViscoTMEM<float>;
template class KITGPI::Modelparameter::ViscoTMEM<double>;
