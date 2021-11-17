#include "ViscoSH.hpp"
#include "../IO/IO.hpp"

using namespace scai;

/*! \brief estimate sum of the memory of all model parameters
 *
 \param dist Distribution
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::ViscoSH<ValueType>::estimateMemory(dmemo::DistributionPtr dist)
{
    /* 6 Parameter in ViscoSH modeling:  rho, Vs, invRho, sWaveModulus, sWaveModulusXZ, sWaveModulusYZ, tauS, tauSXZ, tauSYZ */
    IndexType numParameter = 9;
    return (this->getMemoryUsage(dist, numParameter));
}

/*! \brief Prepare modelparameter for modelling
 *
 * Refreshes the modulus, calculates inverse density and average Values on staggered grid
 *
 \param modelCoordinates Coordinate class object
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters\n");

    // refreshModulus();
    this->getSWaveModulus();
    this->getInverseDensity();
    //sWaveModulus is not needed after averaging. The memory should be freed after calculating the averaging.
    initializeMatrices(dist, ctx, modelCoordinates, comm);
    calculateAveraging();
    purgeMatrices();
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::applyThresholds(Configuration::Configuration const &config)
{
    lama::DenseVector<ValueType> maskS(velocityS); //mask to restore vacuum
    maskS.unaryOp(maskS, common::UnaryOp::SIGN);
    maskS.unaryOp(maskS, common::UnaryOp::ABS);

    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("lowerVSTh"), config.get<ValueType>("lowerVSTh"), 1);
    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("upperVSTh"), config.get<ValueType>("upperVSTh"), 2);
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty

    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("lowerDensityTh"), config.get<ValueType>("lowerDensityTh"), 1);
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("upperDensityTh"), config.get<ValueType>("upperDensityTh"), 2);
    dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed
    dirtyFlagAveraging = true;      // If S-Wave velocity will be changed, averaging needs to be redone

    if (config.getAndCatch("inversionType", 1) == 3 || config.getAndCatch("parameterisation", 0) == 1 || config.getAndCatch("parameterisation", 0) == 2) {
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("lowerPorosityTh", 0.0), config.getAndCatch("lowerPorosityTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(porosity, config.getAndCatch("upperPorosityTh", 1.0), config.getAndCatch("upperPorosityTh", 1.0), 2);

        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("lowerSaturationTh", 0.0), config.getAndCatch("lowerSaturationTh", 0.0), 1);
        Common::searchAndReplace<ValueType>(saturation, config.getAndCatch("upperSaturationTh", 1.0), config.getAndCatch("upperSaturationTh", 1.0), 2);
    }
    velocityS *= maskS;
    density *= maskS;
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
void KITGPI::Modelparameter::ViscoSH<ValueType>::getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    auto distBig = velocityS.getDistributionPtr();
    auto dist = modelPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    
    lama::DenseVector<ValueType> temp;
        
    temp = shrinkMatrix * velocityS;
    modelPerShot.setVelocityS(temp);
    
    temp = shrinkMatrix * tauS;
    modelPerShot.setTauS(temp);
    
    temp = shrinkMatrix * density;
    modelPerShot.setDensity(temp);
    
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
void KITGPI::Modelparameter::ViscoSH<ValueType>::setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    auto distBig = density.getDistributionPtr();
    auto dist = modelPerShot.getDensity().getDistributionPtr();

    scai::lama::CSRSparseMatrix<ValueType> shrinkMatrix = this->getShrinkMatrix(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate);
    shrinkMatrix.assignTranspose(shrinkMatrix);
    
    scai::lama::SparseVector<ValueType> eraseVector = this->getEraseVector(dist, distBig, modelCoordinates, modelCoordinatesBig, cutCoordinate, boundaryWidth);
    scai::lama::SparseVector<ValueType> restoreVector;
    restoreVector = 1.0 - eraseVector;
    
    lama::DenseVector<ValueType> temp;
    
    temp = shrinkMatrix * modelPerShot.getVelocityS(); //transform pershot into big model
    temp *= restoreVector;
    velocityS *= eraseVector;
    velocityS += temp; //take over the values
  
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
KITGPI::Modelparameter::ViscoSH<ValueType>::ViscoSH(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    equationType = "viscosh";
    init(config, ctx, dist, modelCoordinates);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_ASSERT(config.get<IndexType>("ModelRead") != 2, "Read variable model not available for ViscoSH, variable grid is not available here!")
    
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

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Reading model parameter (ViscoSH) from file...\n");

        initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));
        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("FileFormat"));

        HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityS"), config.get<ValueType>("rho"), config.get<ValueType>("tauS"), numRelaxationMechanisms_in, relaxationFrequency_in, config.get<ValueType>("CenterFrequencyCPML"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given values.
 \param ctx Context
 \param dist Distribution
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType>::ViscoSH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    equationType = "viscosh";
    init(ctx, dist, velocityS_const, rho_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given values.
 \param ctx Context
 \param dist Distribution
 \param velocityS_const S-wave velocity given as Scalar
 \param rho_const Density given as Scalar
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in, centerFrequencyCPML_in);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho_const);
    this->initModelparameter(tauS, ctx, dist, tauS_const);
    this->initModelparameter(porosity, ctx, dist, 0.0);
    this->initModelparameter(saturation, ctx, dist, 0.0);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename base filenam of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType>::ViscoSH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    equationType = "viscosh";
    init(ctx, dist, filename, fileFormat);
}

/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename base filename of the model
 \param fileFormat Input file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    this->initModelparameter(velocityS, ctx, dist, filename + ".vs", fileFormat);
    this->initModelparameter(density, ctx, dist, filename + ".density", fileFormat);
    this->initModelparameter(tauS, ctx, dist, filename + ".tauS", fileFormat);
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
KITGPI::Modelparameter::ViscoSH<ValueType>::ViscoSH(const ViscoSH &rhs)
{
    equationType = rhs.equationType;
    sWaveModulus = rhs.sWaveModulus;
    velocityS = rhs.velocityS;
    density = rhs.density;
    inverseDensity = rhs.inverseDensity;
    tauS = rhs.tauS;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
    
    porosity = rhs.porosity;
    saturation = rhs.saturation;
}

/*! \brief Write model to an external file
 *base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::write(std::string filename, scai::IndexType fileFormat) const
{
    IO::writeVector(density, filename + ".density", fileFormat);
    IO::writeVector(velocityS, filename + ".vs", fileFormat);
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
void KITGPI::Modelparameter::ViscoSH<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat)
{
    scai::lama::DenseVector<ValueType> vsRockMatrix;
    
    this->calcVelocityFromModulus(shearModulusRockMatrix, densityRockMatrix, vsRockMatrix);
    
    IO::writeVector(densityRockMatrix, filename + ".densityma", fileFormat);
    IO::writeVector(vsRockMatrix, filename + ".vsma", fileFormat);
};

/*! \brief calculate densityRockMatrix and shearModulusRockMatrix from porosity and saturation */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::calcRockMatrixParameter(Configuration::Configuration const &config)
{        
    scai::lama::DenseVector<ValueType> rho_ma;
    scai::lama::DenseVector<ValueType> mu_ma; 
    scai::lama::DenseVector<ValueType> rho_sat; 
    scai::lama::DenseVector<ValueType> mu_sat;    
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    
    rho_sat = this->getDensity();  
    mu_sat = this->getSWaveModulus();  
    
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
    
    this->setDensityRockMatrix(rho_ma);
    this->setShearModulusRockMatrix(mu_ma);
    
    // initialisation of necessary parameters for gradientCalculation in which a const model is used.
    this->getBiotCoefficient();
}

/*! \brief transform porosity and saturation to Seismic parameters */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::calcWaveModulusFromPetrophysics()
{
    scai::lama::DenseVector<ValueType> densitytemp;
    scai::lama::DenseVector<ValueType> velocityStemp;
    scai::lama::DenseVector<ValueType> mu_sat;
    scai::lama::DenseVector<ValueType> mu_ma;  
    scai::lama::DenseVector<ValueType> beta;
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    
    mu_ma = this->getShearModulusRockMatrix();  

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
    
    this->setDensity(densitytemp);
    this->setVelocityS(velocityStemp);
}

/*! \brief transform Seismic parameters to porosity and saturation  */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::calcPetrophysicsFromWaveModulus()
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
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoSH.initializeMatrices")
        HOST_PRINT(comm, "", "Preparation of the Average Matrix\n");
        // reuse of density average for the swavemodulus : sxz and syz are on the same spot as vx and vy in acoustic modeling
        this->calcAverageMatrixX(modelCoordinates, dist);
        this->calcAverageMatrixY(modelCoordinates, dist);

        averageMatrixX.setContextPtr(ctx);
        averageMatrixY.setContextPtr(ctx);
    }
}

//! \brief Purge Averaging matrices to free memory
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::purgeMatrices()
{
    averageMatrixX.purge();
    averageMatrixY.purge();
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging) {
        SCAI_REGION("Modelparameter.ViscoSH.calculateAveraging")
        // sxz and syz are on the same spot as vx and vy in acoustic modeling
        this->calcAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, averageMatrixX);
        this->calcAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, averageMatrixY);
        this->calcAveragedParameter(tauS, tauSAverageXZ, averageMatrixX);
        this->calcAveragedParameter(tauS, tauSAverageYZ, averageMatrixY);
        dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an viscosh simulation")
    }
    if (relaxationFrequency_in[0] <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an viscosh simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
    centerFrequencyCPML = centerFrequencyCPML_in;
}

/*! \brief calculate reflectivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT)
{
    scai::lama::DenseVector<ValueType> impedance;
    scai::lama::DenseVector<ValueType> impedanceAverageY;
    auto dist = density.getDistributionPtr();
    auto ctx = density.getContextPtr();
    auto const &Dyf = derivatives.getDyf();
    this->calcAverageMatrixY(modelCoordinates, dist);
    averageMatrixY.setContextPtr(ctx);
    impedance = density * velocityS;
    this->calcAveragedParameter(impedance, impedanceAverageY, averageMatrixY);
    averageMatrixY.purge();
    impedanceAverageY *= 2;
    reflectivity = Dyf * impedance;
    reflectivity *= -modelCoordinates.getDH() / DT;
    reflectivity /= impedanceAverageY;
}

/*! \brief Get equationType (viscosh)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::ViscoSH<ValueType>::getEquationType() const
{
    return (equationType);
}

/*! \brief Get reference to velocityP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getVelocityP() const
{
    COMMON_THROWEXCEPTION("There is no velocityP parameter in an sh modelling")
    return (velocityP);
}

/*! \brief Get reference to pWaveModulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getPWaveModulus()
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus parameter in an sh modelling")
    return (pWaveModulus);
}

/*! \brief Get reference to pWaveModulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getPWaveModulus() const
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus parameter in an sh modelling")
    return (pWaveModulus);
}

/*! \brief Get const reference to S-wave modulus (viscoelastic case)
 *
 * if S-Wave Modulus is dirty (eg. because of changes in velocityS, the S-Wave modulus will be recalculated
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getSWaveModulus()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "S-Wave Modulus has to be recalculated! ");
    return (sWaveModulus);
}

/*! \brief Get reference to bulkModulusRockMatrix
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getBulkModulusRockMatrix() const
{
    COMMON_THROWEXCEPTION("There is no bulkModulusRockMatrix parameter in an sh modelling")
    return (bulkModulusRockMatrix);
}

/*! \brief Get reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tau parameter in an sh modelling")
    return (tauP);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXY);
}

/*! \brief Get reference to tauS xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ViscoSH<ValueType>::getTauSAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an sh modelling")
    return (tauSAverageXY);
}

/*! \brief Overloading * Operation
 *
 \param rhs ValueType factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> KITGPI::Modelparameter::ViscoSH<ValueType>::operator*(ValueType rhs)
{
    KITGPI::Modelparameter::ViscoSH<ValueType> result(*this);
    result *= rhs;
    return result;
}

/*! \brief free function to multiply
 *
 \param lhs ValueType factor with which the vectors are multiplied.
 \param rhs Vector
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    return rhs * lhs;
}

/*! \brief Overloading *= Operation
 *
 \param rhs valueType factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> &KITGPI::Modelparameter::ViscoSH<ValueType>::operator*=(ValueType const &rhs)
{
    density *= rhs;
    velocityS *= rhs;
    tauS *= rhs;
    porosity *= rhs;
    saturation *= rhs;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> KITGPI::Modelparameter::ViscoSH<ValueType>::operator+(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoSH<ValueType> result(*this);
    result += rhs;
    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> &KITGPI::Modelparameter::ViscoSH<ValueType>::operator+=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    density += rhs.density;
    velocityS += rhs.velocityS;
    tauS += rhs.tauS;
    porosity += rhs.porosity;
    saturation += rhs.saturation;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> KITGPI::Modelparameter::ViscoSH<ValueType>::operator-(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    KITGPI::Modelparameter::ViscoSH<ValueType> result(*this);
    result -= rhs;
    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> &KITGPI::Modelparameter::ViscoSH<ValueType>::operator-=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    density -= rhs.density;
    velocityS -= rhs.velocityS;
    tauS -= rhs.tauS;
    porosity -= rhs.porosity;
    saturation -= rhs.saturation;

    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ViscoSH<ValueType> &KITGPI::Modelparameter::ViscoSH<ValueType>::operator=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs)
{
    velocityS = rhs.velocityS;
    tauS = rhs.tauS;
    density = rhs.density;
    porosity = rhs.porosity;
    saturation = rhs.saturation;
    
    shearModulusRockMatrix = rhs.shearModulusRockMatrix;
    densityRockMatrix = rhs.densityRockMatrix;
    
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityS = rhs.getVelocityS();
    tauS = rhs.getTauS();
    density = rhs.getDensity();
    porosity = rhs.getPorosity();
    saturation = rhs.getSaturation();
    
    shearModulusRockMatrix = rhs.getShearModulusRockMatrix();
    densityRockMatrix = rhs.getDensityRockMatrix();
    
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading -= Operation (called in base class)
 *
 \param rhs Abstract model which is substracted.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityS -= rhs.getVelocityS();
    tauS -= rhs.getTauS();
    density -= rhs.getDensity();
    porosity -= rhs.getPorosity();
    saturation -= rhs.getSaturation();
    
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

/*! \brief function for overloading += Operation (called in base class)
 *
 \param rhs Abstract model which is added.
 */
template <typename ValueType>
void KITGPI::Modelparameter::ViscoSH<ValueType>::plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    velocityS += rhs.getVelocityS();
    tauS += rhs.getTauS();
    density += rhs.getDensity();
    porosity += rhs.getPorosity();
    saturation += rhs.getSaturation();
    
    dirtyFlagInverseDensity = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::ViscoSH<float>;
template class KITGPI::Modelparameter::ViscoSH<double>;
