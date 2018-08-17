#include "Viscoelastic.hpp"
using namespace scai;



/*! \brief Prepare modellparameter for visco-elastic modelling
 *
 * Applies Equation 12 from Bohlen 2002 and refreshes the modulus
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "Preparation of the model parameters\n"); 
    
    //refreshModulus 
    this->getPWaveModulus();  
    this->getSWaveModulus();
    initializeMatrices(dist, ctx, config, comm);
  
    this->getInverseDensity();
    calculateAveraging();
    

    HOST_PRINT(comm, "Model ready!\n\n");
}

/*! \brief Apply thresholds to model parameters
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::applyThresholds(Configuration::Configuration const &config) {
    lama::DenseVector<ValueType> maskP(velocityP); //mask to restore vacuum
    maskP.unaryOp(maskP,common::UnaryOp::SIGN);
    maskP.unaryOp(maskP,common::UnaryOp::ABS);
    
    lama::DenseVector<ValueType> maskS(velocityS); //mask to restore acoustic media
    maskS.unaryOp(maskP,common::UnaryOp::SIGN);
    maskS.unaryOp(maskP,common::UnaryOp::ABS);
    
    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("lowerVPTh"), config.get<ValueType>("lowerVPTh"), 1);
    Common::searchAndReplace<ValueType>(velocityP, config.get<ValueType>("upperVPTh"), config.get<ValueType>("upperVPTh"), 2);
    
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("lowerDensityTh"), config.get<ValueType>("lowerDensityTh"), 1);
    Common::searchAndReplace<ValueType>(density, config.get<ValueType>("upperDensityTh"), config.get<ValueType>("upperDensityTh"), 2);
    
    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("lowerVSTh"), config.get<ValueType>("lowerVSTh"), 1);
    Common::searchAndReplace<ValueType>(velocityS, config.get<ValueType>("upperVSTh"), config.get<ValueType>("upperVSTh"), 2);
     
    velocityP *= maskP;
    density *= maskP;
    velocityS *= maskS;
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    equationType = "viscoelastic";
    init(config, ctx, dist);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    if (config.get<IndexType>("ModelRead")) {

        HOST_PRINT(dist->getCommunicatorPtr(), "Reading model parameter from file...\n");

        init(ctx, dist, config.get<std::string>("ModelFilename"), config.get<IndexType>("PartitionedIn"));
        initRelaxationMechanisms(config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));

        HOST_PRINT(dist->getCommunicatorPtr(), "Finished with reading of the model parameter!\n\n");

    } else {
        init(ctx, dist, config.get<ValueType>("velocityP"), config.get<ValueType>("velocityS"), config.get<ValueType>("rho"), config.get<ValueType>("tauP"), config.get<ValueType>("tauS"), config.get<IndexType>("numRelaxationMechanisms"), config.get<ValueType>("relaxationFrequency"));
    }

    if (config.get<IndexType>("ModelWrite")) {
        write(config.get<std::string>("ModelFilename") + ".out", config.get<IndexType>("PartitionedOut"));
	std::cout << "been here\n\n";
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
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    equationType = "viscoelastic";
    init(ctx, dist, velocityP_const, velocityS_const, rho_const, tauP_const, tauS_const, numRelaxationMechanisms_in, relaxationFrequency_in);
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
}

/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the five given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 \param tauP_const TauP given as Scalar
 \param tauS_const TauS given as Scalar
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityP_const, ValueType velocityS_const, ValueType rho_const, ValueType tauP_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    this->initModelparameter(velocityP, ctx, dist, velocityP_const);
    this->initModelparameter(velocityS, ctx, dist, velocityS_const);
    this->initModelparameter(density, ctx, dist, rho_const);
    this->initModelparameter(tauS, ctx, dist, tauS_const);
    this->initModelparameter(tauP, ctx, dist, tauP_const);
}

/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    equationType = "viscoelastic";
    init(ctx, dist, filename, partitionedIn);
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedIn Partitioned input
 *
 *  Calculates pWaveModulus with
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    std::string filenameVelocityP = filename + ".vp.mtx";
    std::string filenameVelocityS = filename + ".vs.mtx";
    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";

    this->initModelparameter(velocityS, ctx, dist, filenameVelocityS, partitionedIn);
    this->initModelparameter(velocityP, ctx, dist, filenameVelocityP, partitionedIn);
    this->initModelparameter(density, ctx, dist, filenamedensity, partitionedIn);
    this->initModelparameter(tauS, ctx, dist, filenameTauS, partitionedIn);
    this->initModelparameter(tauP, ctx, dist, filenameTauP, partitionedIn);
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
}


/*! \brief Write model to an external file
 *
 \param filename For the P-wave velocity ".vp.mtx" is added, for the S-wave velocity ".vs.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{

    std::string filenamedensity = filename + ".density.mtx";
    std::string filenameTauP = filename + ".tauP.mtx";
    std::string filenameTauS = filename + ".tauS.mtx";
    std::string filenameP = filename + ".vp.mtx";
    std::string filenameS = filename + ".vs.mtx";
    
    this->writeModelparameter(density, filenamedensity, partitionedOut);
    this->writeModelparameter(tauP, filenameTauP, partitionedOut);
    this->writeModelparameter(tauS, filenameTauS, partitionedOut);
    this->writeModelparameter(velocityP, filenameP, partitionedOut);
    this->writeModelparameter(velocityS, filenameS, partitionedOut);

};

//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration config, scai::dmemo::CommunicatorPtr comm)
{
    initializeMatrices(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"), config.get<ValueType>("DT"), comm);
}

//! \brief Initializsation of the Averaging matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, scai::dmemo::CommunicatorPtr /*comm*/)
{
    if (dirtyFlagAveraging)
    {
    SCAI_REGION("initializeMatrices")

    this->calcDensityAverageMatrixX(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixY(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixZ(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixXY(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixXZ(NX, NY, NZ, dist);
    this->calcSWaveModulusAverageMatrixYZ(NX, NY, NZ, dist);
    
    DensityAverageMatrixX.setContextPtr(ctx);
    DensityAverageMatrixY.setContextPtr(ctx);
    DensityAverageMatrixZ.setContextPtr(ctx);
    sWaveModulusAverageMatrixXY.setContextPtr(ctx);
    sWaveModulusAverageMatrixXZ.setContextPtr(ctx);
    sWaveModulusAverageMatrixYZ.setContextPtr(ctx);
    }
}

/*! \brief calculate averaged vectors
 *
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::calculateAveraging()
{
    if (dirtyFlagAveraging)
    {
    this->calculateInverseAveragedDensity(density, inverseDensityAverageX, DensityAverageMatrixX);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageY, DensityAverageMatrixY);
    this->calculateInverseAveragedDensity(density, inverseDensityAverageZ, DensityAverageMatrixZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedSWaveModulus(sWaveModulus, sWaveModulusAverageYZ, sWaveModulusAverageMatrixYZ);
    this->calculateAveragedTauS(tauS, tauSAverageXY, sWaveModulusAverageMatrixXY);
    this->calculateAveragedTauS(tauS, tauSAverageXZ, sWaveModulusAverageMatrixXZ);
    this->calculateAveragedTauS(tauS, tauSAverageYZ, sWaveModulusAverageMatrixYZ);
    dirtyFlagAveraging = false;
    }
}

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    if (numRelaxationMechanisms_in < 1) {
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-elastic simulation")
    }
    if (relaxationFrequency_in <= 0) {
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-elastic simulation")
    }
    numRelaxationMechanisms = numRelaxationMechanisms_in;
    relaxationFrequency = relaxationFrequency_in;
}

/*! \brief Get equationType (viscoelastic)
 */
template <typename ValueType>
std::string KITGPI::Modelparameter::Viscoelastic<ValueType>::getEquationType() const
{
    return (equationType);
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
	HOST_PRINT(velocityP.getDistributionPtr()->getCommunicatorPtr(), "P-Wave modulus will be calculated from density,velocityP,tauP and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityP,density,pWaveModulus);
	/* Set circular frequency w = 2 * pi * relaxation frequency */
	ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
	ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency);
	
	ValueType sum = w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);

	/* Scaling the P-wave Modulus */

	auto temp = lama::eval<lama::DenseVector<ValueType>>( 1.0 + sum * tauP );
	pWaveModulus = pWaveModulus / temp;
    dirtyFlagPWaveModulus = false;
    }
    
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
	HOST_PRINT(velocityS.getDistributionPtr()->getCommunicatorPtr(), "S-Wave modulus will be calculated from density,velocityS,tauS and relaxationFrequency \n");
        this->calcModulusFromVelocity(velocityS,density,sWaveModulus);
	/* Set circular frequency w = 2 * pi * relaxation frequency */
	ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
	ValueType tauSigma = 1.0 / (2.0 * M_PI * relaxationFrequency);
	
	ValueType sum = w_ref * w_ref * tauSigma * tauSigma / (1.0 + w_ref * w_ref * tauSigma * tauSigma);

	/* Scaling the S-wave Modulus */
	auto temp = lama::eval<lama::DenseVector<ValueType>>( 1.0 + sum * tauS );
	sWaveModulus = sWaveModulus / temp;
        dirtyFlagSWaveModulus = false;
    }
    
    return (pWaveModulus);
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
    
    dirtyFlagInverseDensity=true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging=true;
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
	
    dirtyFlagInverseDensity=true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging=true;
    
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
    density = density -= rhs.density;
    tauS -= rhs.tauS;
    tauP -= rhs.tauP;
    velocityP -= rhs.velocityP;
    velocityS -= rhs.velocityS;
    
    dirtyFlagInverseDensity=true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging=true;
    return *this;
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> &KITGPI::Modelparameter::Viscoelastic<ValueType>::operator=(KITGPI::Modelparameter::Viscoelastic<ValueType> const &rhs)
{
    pWaveModulus = rhs.pWaveModulus;
    sWaveModulus = rhs.sWaveModulus;
    velocityP = rhs.velocityP;
    velocityS = rhs.velocityS;
    density = rhs.density;
    inverseDensity = rhs.inverseDensity;
    tauS = rhs.tauS;
    tauP = rhs.tauP;
    relaxationFrequency = rhs.relaxationFrequency;
    numRelaxationMechanisms = rhs.numRelaxationMechanisms;
    dirtyFlagInverseDensity = rhs.dirtyFlagInverseDensity;    
    dirtyFlagPWaveModulus = rhs.dirtyFlagPWaveModulus;
    dirtyFlagSWaveModulus = rhs.dirtyFlagSWaveModulus;
    dirtyFlagAveraging= rhs.dirtyFlagAveraging;
    return *this;
}

/*! \brief function for overloading = Operation (called in base class)
 *
 \param rhs Abstract model which is assigned.
 */
template <typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    pWaveModulus = rhs.getPWaveModulus();
    sWaveModulus = rhs.getSWaveModulus();
    velocityP = rhs.getVelocityP();
    velocityS = rhs.getVelocityS();
    inverseDensity = rhs.getInverseDensity();
    density = rhs.getDensity();
    tauS = rhs.getTauS();
    tauP = rhs.getTauP();
    relaxationFrequency = rhs.getRelaxationFrequency();
    numRelaxationMechanisms = rhs.getNumRelaxationMechanisms();
    dirtyFlagInverseDensity = rhs.getDirtyFlagInverseDensity();
    dirtyFlagPWaveModulus = rhs.getDirtyFlagPWaveModulus();
    dirtyFlagSWaveModulus= rhs.getDirtyFlagSWaveModulus();
    dirtyFlagAveraging = rhs.getDirtyFlagAveraging();
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
    dirtyFlagInverseDensity = true;
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagAveraging = true;
}

template class KITGPI::Modelparameter::Viscoelastic<float>;
template class KITGPI::Modelparameter::Viscoelastic<double>;
