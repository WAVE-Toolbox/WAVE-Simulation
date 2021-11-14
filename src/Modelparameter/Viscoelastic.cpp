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
    /* 15 Parameter in Viscoelastic modeling:  rho, Vp, Vs, invRhoX,invRhoY, invRhoZ,  bulk modulus, sWaveModulus, sWaveModulusXY, sWaveModulusXZ, sWaveModulusYZ, tauP, tauS, tauSXY, tauSXZ, tauSYZ*/
    IndexType numParameter = 15;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modelparameter for viscoelastic modelling
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

    if (config.getAndCatch("inversionType", 1) == 3 || config.getAndCatch("parameterisation", 0) == 1 || config.getAndCatch("parameterisation", 0) == 2) {
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = density.getDistributionPtr();
    auto dist = modelPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    
    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * velocityP;
    modelPerShot.setVelocityP(temp);
        
    temp = shrinkMatrix * velocityS;
    modelPerShot.setVelocityS(temp);
    
    temp = shrinkMatrix * density;
    modelPerShot.setDensity(temp);
    
    temp = shrinkMatrix * tauP;
    modelPerShot.setTauP(temp);
        
    temp = shrinkMatrix * tauS;
    modelPerShot.setTauS(temp);
    
    temp = shrinkMatrix * porosity;
    modelPerShot.setPorosity(temp);
    
    temp = shrinkMatrix * saturation;
    modelPerShot.setSaturation(temp);
}

/*! \brief If stream configuration is used, get a pershot model from the big model
 \param modelPerShot pershot model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate cut coordinate 
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
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

    temp = shrinkMatrix * modelPerShot.getTauP(); //transform pershot into big model
    temp *= restoreVector;
    tauP *= eraseVector;
    tauP += temp; //take over the values
  
    temp = shrinkMatrix * modelPerShot.getTauS(); //transform pershot into big model
    temp *= restoreVector;
    tauS *= eraseVector;
    tauS += temp; //take over the values
    
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
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2, "Read variable model not available for Viscoelastic, variable grid is not available here!")

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

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model (viscoelastic) parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));
        initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"), config.get<ValueType>("tauP"), config.get<ValueType>("tauS"), numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));
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
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    equationType = "viscoelastic";
    init(ctx, dist, velocityP_const, velocityS_const, rho_const, tauP_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the five given scalar values.
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho_const);
    this->initModelparameter(tauS, ctx, dist, tauS_const);
    this->initModelparameter(tauP, ctx, dist, tauP_const);
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
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
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
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Write model to an external file
 *base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    IO::writeVector(density, filename + ".density", fileFormat);
    IO::writeVector(velocityP, filename + ".vp", fileFormat);
    IO::writeVector(velocityS, filename + ".vs", fileFormat);
    IO::writeVector(tauP, filename + ".tauP", fileFormat);
    IO::writeVector(tauS, filename + ".tauS", fileFormat);
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat)
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

/*! \brief calculate densityRockMatrix and bulkModulusRockMatrix from porosity and saturation */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calcRockMatrixParameter(Configuration::Configuration const &config)
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
     
    rho_ma = scai::lama::cast<ValueType>(rho_ma);
    Common::replaceInvalid<ValueType>(rho_ma, 0.0);
    
    // calculate mu_ma
    temp1 = 1 - beta;
    temp1 = 1 / temp1;
    mu_ma = mu_sat * temp1;
     
    mu_ma = scai::lama::cast<ValueType>(mu_ma);
    Common::replaceInvalid<ValueType>(mu_ma, 0.0);
    
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
    K_ma = scai::lama::sqrt(temp1);
    K_ma -= b;
    K_ma /= a;
    K_ma /= 2;
    
    K_ma = scai::lama::cast<ValueType>(K_ma);    
    Common::replaceInvalid<ValueType>(K_ma, 0.0);
    
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calcWaveModulusFromPetrophysics()
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
     
    // calculate density
    temp1 = saturation * DensityWater;
    temp2 = 1 - saturation;
    temp2 *= DensityAir;
    temp1 += temp2;
    temp1 *= porosity; 
    temp2 = 1 - porosity;
    densitytemp = temp2 * this->getDensityRockMatrix();
    densitytemp += temp1;
    
    densitytemp = scai::lama::cast<ValueType>(densitytemp);    
    Common::replaceInvalid<ValueType>(densitytemp, 0.0);
  
    // calculate vs
    temp1 = 1 - beta;    
    mu_sat = mu_ma * temp1;    
    velocityStemp = mu_sat / densitytemp;
    velocityStemp = scai::lama::sqrt(velocityStemp);
     
    velocityStemp = scai::lama::cast<ValueType>(velocityStemp);
    Common::replaceInvalid<ValueType>(velocityStemp, 0.0);
        
    // calculate vp
    M = this->getBulkModulusM();
    temp2 = M * beta;
    temp2 *= beta;
    temp1 = 1 - beta;
    K_sat = K_ma * temp1;
    K_sat += temp2;
    velocityPtemp = 4 / 3 * mu_sat;
    velocityPtemp += K_sat;    
    velocityPtemp /= densitytemp;
    velocityPtemp = scai::lama::sqrt(velocityPtemp);  
    
    velocityPtemp = scai::lama::cast<ValueType>(velocityPtemp);
    Common::replaceInvalid<ValueType>(velocityPtemp, 0.0);  
    
    this->setDensity(densitytemp);
    this->setVelocityS(velocityStemp);
    this->setVelocityP(velocityPtemp);
}

/*! \brief transform Seismic parameters to porosity and saturation  */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calcPetrophysicsFromWaveModulus()
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
    
    porositytemp = scai::lama::cast<ValueType>(porositytemp); 
    Common::replaceInvalid<ValueType>(porositytemp, 0.0);
    
    this->setPorosity(porositytemp);
}

//! \brief Initialization of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates coordinate class object
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.Viscoelastic.initializeMatrices")
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
        SCAI_REGION("Modelparameter.Viscoelastic.calculateAveraging")
        this->calcInverseAveragedParameter(density, inverseDensityAverageX, averageMatrixX);
        this->calcInverseAveragedParameter(density, inverseDensityAverageY, averageMatrixY);
        this->calcInverseAveragedParameter(density, inverseDensityAverageZ, averageMatrixZ);
        this->calcAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, averageMatrixXY);
        this->calcAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, averageMatrixXZ);
        this->calcAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, averageMatrixYZ);
        this->calcAveragedParameter(tauS, tauSAverageXY, averageMatrixXY);
        this->calcAveragedParameter(tauS, tauSAverageXZ, averageMatrixXZ);
        this->calcAveragedParameter(tauS, tauSAverageYZ, averageMatrixYZ);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an viscoelastic simulation")
    }
    if (relaxationFrequency_in[0] <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an viscoelastic simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
    centerFrequencyCPML = centerFrequencyCPML_in;
}

/*! \brief calculate reflectivity from permittivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT)
{
    scai::lama::DenseVector<ValueType> impedance;
    scai::lama::DenseVector<ValueType> impedanceAverageY;
    auto dist = density.getDistributionPtr();
    auto ctx = density.getContextPtr();
    auto const &Dyf = derivatives.getDyf();
    this->calcAverageMatrixY(modelCoordinates, dist);
    averageMatrixY.setContextPtr(ctx);
    impedance = density * velocityP;
    this->calcAveragedParameter(impedance, impedanceAverageY, averageMatrixY);
    averageMatrixY.purge();
    impedanceAverageY *= 2;
    reflectivity = Dyf * impedance;
    reflectivity *= -modelCoordinates.getDH() / DT;
    reflectivity /= impedanceAverageY;
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
    COMMON_THROWEXCEPTION("Inverse density is not set for viscoelastic modelling")
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
        // The input vp is defined at the reference frequency, so the P-wave modulus must be scaled to the relaxed value where frequency equals zero, see eq.12 in Bohlen, 2002 or eq.2 in Fabien-Ouellet, et al., 2017.
        HOST_PRINT(velocityP.getDistributionPtr()->getCommunicatorPtr(), "P-Wave modulus will be calculated from density, velocityP, tauP, centerFrequencyCPML and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityP, density, pWaveModulus);
        /* Set circular frequency w = 2 * pi * relaxation frequency */
        ValueType w_ref = 2.0 * M_PI * centerFrequencyCPML;
        ValueType sum = 0;
        for (int l=0; l<numRelaxationMechanisms; l++) {
            ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency[l]);

            sum += w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);
        }
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
        // The input vs is defined at the reference frequency, so the S-wave modulus must be scaled to the relaxed value where frequency equals zero, see eq.12 in Bohlen, 2002 or eq.2 in Fabien-Ouellet, et al., 2017.
        HOST_PRINT(velocityS.getDistributionPtr()->getCommunicatorPtr(), "S-Wave modulus will be calculated from density, velocityS, tauS, centerFrequencyCPML and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityS, density, sWaveModulus);
        /* Set circular frequency w = 2 * pi * relaxation frequency */
        ValueType w_ref = 2.0 * M_PI * centerFrequencyCPML;
        ValueType sum = 0;
        for (int l=0; l<numRelaxationMechanisms; l++) {
            ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency[l]);

            sum += w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);
        }
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
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "S-Wave Modulus has to be recalculated! ");
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
    density -= rhs.density;
    tauS -= rhs.tauS;
    tauP -= rhs.tauP;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator=(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    bulkModulusRockMatrix = rhs.bulkModulusRockMatrix;
    shearModulusRockMatrix = rhs.shearModulusRockMatrix;
    densityRockMatrix = rhs.densityRockMatrix;
    
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
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    
    bulkModulusRockMatrix = rhs.getBulkModulusRockMatrix();
    shearModulusRockMatrix = rhs.getShearModulusRockMatrix();
    densityRockMatrix = rhs.getDensityRockMatrix();
    
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityP += rhs.getVelocityP();
    velocityS += rhs.getVelocityS();
    density += rhs.getDensity();
    tauS += rhs.getTauS();
    tauP += rhs.getTauP();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Viscoelastic<float>;
template class KITGPI::Modelparameter::Viscoelastic<double>;
