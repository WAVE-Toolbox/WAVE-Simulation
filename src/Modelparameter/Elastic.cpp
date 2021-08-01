#include "Elastic.hpp"
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

/*! \brief Prepare modelparameter for modelling
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
    maskS.unaryOp(maskS, common::UnaryOp::SIGN);
    maskS.unaryOp(maskS, common::UnaryOp::ABS);

    Common::replaceVpwithVs<ValueType>(velocityP, velocityS, config.get<ValueType>("lowerVpVsRatioTh"), 1);
    Common::replaceVpwithVs<ValueType>(velocityP, velocityS, config.get<ValueType>("upperVpVsRatioTh"), 2);

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

    if (config.get<IndexType>("inversionType") == 3 || config.get<IndexType>("parameterisation") == 1 || config.get<IndexType>("parameterisation") == 2) {
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("lowerPorosityTh", 0.0), config.getAndCatch("lowerPorosityTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("upperPorosityTh", 1.0), config.getAndCatch("upperPorosityTh", 1.0), 2);

        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("lowerSaturationTh", 0.0), config.getAndCatch("lowerSaturationTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("upperSaturationTh", 1.0), config.getAndCatch("upperSaturationTh", 1.0), 2);
    }
    velocityP *= maskP;
    density *= maskP;
    velocityS *= maskS;
    porosity *= maskS;
    saturation *= maskS;
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = density.getDistributionPtr();
    auto dist = modelPerShot.getDensity().getDistributionPtr();
//     auto comm = dist.getCommunicatorPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);

    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * velocityP;
    modelPerShot.setVelocityP(temp);
        
    temp = shrinkMatrix * velocityS;
    modelPerShot.setVelocityS(temp);
    
    temp = shrinkMatrix * density;
    modelPerShot.setDensity(temp);
    
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
void KITGPI::Modelparameter::Elastic<ValueType>::setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    auto distBig = density.getDistributionPtr();
    auto dist = modelPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate, boundaryWidth);
    scai::lama::SparseVector<ValueType> restoreVector;
    restoreVector = 1.0 - eraseVector;
    
    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * modelPerShot.getVelocityP(); //transform pershot into big model
    temp *= restoreVector;
    velocityP *= eraseVector;
    velocityP += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getVelocityS(); //transform pershot into big model
    temp *= restoreVector;
    velocityS *= eraseVector;
    velocityS += temp; //take over the values

    temp = shrinkMatrix * modelPerShot.getDensity(); //transform pershot into big model
    temp *= restoreVector;
    density *= eraseVector;
    density += temp; //take over the values
    
    temp = shrinkMatrix * modelPerShot.getPorosity(); //transform pershot into big model
    temp *= restoreVector;
    porosity *= eraseVector;
    porosity += temp; //take over the values
    
    temp = shrinkMatrix * modelPerShot.getSaturation(); //transform pershot into big model
    temp *= restoreVector;
    saturation *= eraseVector;
    saturation += temp; //take over the values
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
 \param velocityP_const P-wave velocity given as Scalar
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::Elastic<ValueType>::Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const)
{
    equationType = "elastic";
    init(ctx, dist, velocityP_const, velocityS_const, rho_const);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param velocityP_const P-wave velocity given as Scalar
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const)
{
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho_const);
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
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        this->initModelparameter(porosity, ctx, dist, filename + ".porosity", fileFormat);
        this->initModelparameter(saturation, ctx, dist, filename + ".saturation", fileFormat);
    } else {
        this->initModelparameter(porosity, ctx, dist, 0.0);
        this->initModelparameter(saturation, ctx, dist, 0.0);
    }
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
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Write model to an external file
 *
 \param filename base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    IO::writeVector(density, filename + ".density", fileFormat);
    IO::writeVector(velocityP, filename + ".vp", fileFormat);
    IO::writeVector(velocityS, filename + ".vs", fileFormat);
    if (this->getInversionType() == 3 || this->getParameterisation() == 1 || this->getParameterisation() == 2) {
        IO::writeVector(porosity, filename + ".porosity", fileFormat);
        IO::writeVector(saturation, filename + ".saturation", fileFormat);
    }
};

/*! \brief Write model to an external file
 *base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat)
{
    scai::lama::DenseVector<ValueType> vsRockMatrix;
    scai::lama::DenseVector<ValueType> vpRockMatrix;
    scai::lama::DenseVector<ValueType> pWaveModulusRockMatrix;
    
    pWaveModulusRockMatrix = bulkModulusRockMatrix + 4 / 3 * shearModulusRockMatrix;
    this->calcVelocityFromModulus(shearModulusRockMatrix, densityRockMatrix, vsRockMatrix);
    this->calcVelocityFromModulus(pWaveModulusRockMatrix, densityRockMatrix, vpRockMatrix);
    
    IO::writeVector(densityRockMatrix, filename + ".densityma", fileFormat);
    IO::writeVector(vsRockMatrix, filename + ".vsma", fileFormat);
    IO::writeVector(vpRockMatrix, filename + ".vpma", fileFormat);
};

/*! \brief calculate densityRockMatrix, shearModulusRockMatrix and bulkModulusRockMatrix from porosity and saturation */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::calcRockMatrixParameter(Configuration::Configuration const &config)
{    
    scai::lama::DenseVector<ValueType> rho_ma;
    scai::lama::DenseVector<ValueType> mu_ma;  
    scai::lama::DenseVector<ValueType> K_ma;  
    scai::lama::DenseVector<ValueType> Kf;  
    scai::lama::DenseVector<ValueType> a; 
    scai::lama::DenseVector<ValueType> b; 
    scai::lama::DenseVector<ValueType> c; 
    scai::lama::DenseVector<ValueType> rho_sat;
    scai::lama::DenseVector<ValueType> K_sat; 
    scai::lama::DenseVector<ValueType> mu_sat;   
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    ValueType tempValue;
    
    rho_sat = this->getDensity();  
    mu_sat = this->getSWaveModulus();  
    K_sat = this->getPWaveModulus();  
    K_sat -= 4 / 3 * mu_sat;
    Common::searchAndReplace<ValueType>(K_sat, 0.0, 0.0, 1);
    
    // Based on Gassmann equation 
    beta = this->getBiotCoefficient();
    
    // calculate density_ma
    temp1 = saturation * DensityWater;
    temp2 = 1 - saturation;
    temp2 *= DensityAir;
    temp1 += temp2;
    temp1 *= porosity;  
    rho_ma = rho_sat - temp1;
    temp1 = 1 - porosity;
    temp1 = 1 / temp1;
    rho_ma *= temp1;
    
    Common::replaceInvalid<ValueType>(rho_ma, 0.0);
    Common::searchAndReplace<ValueType>(rho_ma, 0.0, 0.0, 1);
    
    // calculate mu_ma
    temp1 = 1 - beta;
    temp1 = 1 / temp1;
    mu_ma = mu_sat * temp1;
     
    Common::replaceInvalid<ValueType>(mu_ma, 0.0);
    Common::searchAndReplace<ValueType>(mu_ma, 0.0, 0.0, 1);
    
    // calculate K_ma
    Kf = this->getBulkModulusKf();
    temp1 = 1 - beta;
    a = temp1 / Kf;
    tempValue = 1 / CriticalPorosity - 1;
    c = -tempValue * K_sat;
    b = temp1 * tempValue;
    temp1 = beta / CriticalPorosity;
    b += temp1;
    temp1 = K_sat / Kf;
    b -= temp1;
    temp1 = b * b;
    temp2 = 4 * a * c;
    temp1 -= temp2;    
    Common::searchAndReplace<ValueType>(temp1, 0.0, 0.0, 1);
    K_ma = scai::lama::sqrt(temp1);
    K_ma -= b;
    K_ma /= a;
    K_ma /= 2;
       
    Common::replaceInvalid<ValueType>(K_ma, 0.0);
    Common::searchAndReplace<ValueType>(K_ma, 0.0, 0.0, 1);
    
    this->setDensityRockMatrix(rho_ma);
    this->setShearModulusRockMatrix(mu_ma);
    this->setBulkModulusRockMatrix(K_ma);
    
    // initialisation of necessary parameters for gradientCalculation in which a const model is used.
    this->getBiotCoefficient();
    this->getBulkModulusKf();
    this->getBulkModulusM();
}

/*! \brief transform porosity and saturation to Seismic parameters */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::calcWaveModulusFromPetrophysics()
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> velocityPtemp;
    scai::lama::DenseVector<ValueType> K_sat;
    scai::lama::DenseVector<ValueType> mu_sat;
    scai::lama::DenseVector<ValueType> K_ma;
    scai::lama::DenseVector<ValueType> mu_ma;
    scai::lama::DenseVector<ValueType> M;    
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    
    mu_ma = this->getShearModulusRockMatrix();  
    K_ma = this->getBulkModulusRockMatrix();  

    // Based on Gassmann equation
    beta = this->getBiotCoefficient();
     
//     std::cout << "density\n" << std::endl;
    // calculate density
    temp1 = saturation * DensityWater;
    temp2 = 1 - saturation;
    temp2 *= DensityAir;
    temp1 += temp2;
    temp1 *= porosity; 
    temp2 = 1 - porosity;
    densitytemp = temp2 * this->getDensityRockMatrix();
    densitytemp += temp1;
     
    Common::replaceInvalid<ValueType>(densitytemp, 0.0);
    Common::searchAndReplace<ValueType>(densitytemp, 0.0, 0.0, 1);
  
    // calculate vs
//     std::cout << "vs\n" << std::endl;
    temp1 = 1 - beta;    
    mu_sat = mu_ma * temp1;    
    Common::searchAndReplace<ValueType>(mu_sat, 0.0, 0.0, 1);
    velocityStemp = mu_sat / densitytemp;
    velocityStemp = scai::lama::sqrt(velocityStemp);
     
    Common::replaceInvalid<ValueType>(velocityStemp, 0.0);
        
//     std::cout << "vp\n" << std::endl;
    // calculate vp
    M = this->getBulkModulusM();
    temp2 = M * beta;
    temp2 *= beta;
    temp1 = 1 - beta;
    K_sat = K_ma * temp1;
    K_sat += temp2;
    Common::searchAndReplace<ValueType>(K_sat, 0.0, 0.0, 1);
    velocityPtemp = 4 / 3 * mu_sat;
    velocityPtemp += K_sat;    
    velocityPtemp /= densitytemp;
    velocityPtemp = scai::lama::sqrt(velocityPtemp);  
    
    Common::replaceInvalid<ValueType>(velocityPtemp, 0.0);
    
//     std::cout << "finished\n" << std::endl;
    this->setDensity(densitytemp);
    this->setVelocityS(velocityStemp);
    this->setVelocityP(velocityPtemp);
}

/*! \brief transform Seismic parameters to porosity and saturation  */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::calcPetrophysicsFromWaveModulus()
{
    scai::lama::DenseVector<ValueType> mu_sat;
    scai::lama::DenseVector<ValueType> mu_ma;
    scai::lama::DenseVector<ValueType> porositytemp;
    
    mu_sat = this->getSWaveModulus();  
    mu_ma = this->getShearModulusRockMatrix();  

    // Based on Gassmann equation
    porositytemp = mu_sat / mu_ma;
    porositytemp = 1 - porositytemp;
    porositytemp *= CriticalPorosity;
    
//     porositytemp = scai::lama::cast<ValueType>(porositytemp); 
    Common::replaceInvalid<ValueType>(porositytemp, 0.0);
    
    this->setPorosity(porositytemp);
}

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Elastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.Elastic.initializeMatrices")
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
    porosity *= rhs;
    saturation *= rhs;

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
    porosity += rhs.porosity;
    saturation += rhs.saturation;

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
    density -= rhs.density;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

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
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    bulkModulusRockMatrix = rhs.bulkModulusRockMatrix;
    shearModulusRockMatrix = rhs.shearModulusRockMatrix;
    densityRockMatrix = rhs.densityRockMatrix;
    std::cout << "operator=\n" << std::endl;
    
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
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    
    bulkModulusRockMatrix = rhs.getBulkModulusRockMatrix();
    shearModulusRockMatrix = rhs.getShearModulusRockMatrix();
    densityRockMatrix = rhs.getDensityRockMatrix();
//     std::cout << "assign\n" << std::endl;
    
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
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    
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
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Elastic<float>;
template class KITGPI::Modelparameter::Elastic<double>;
