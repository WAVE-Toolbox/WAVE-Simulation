
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "Modelparameter.hpp"
#include "../PartitionedInOut/PartitionedInOut.hpp"

namespace KITGPI {
    
    //! \brief Modelparameter namespace
    namespace Modelparameter {
        
        //! Class for Modelparameter for acoustic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the acoustic finite-difference simulation.
         */
        template<typename ValueType>
        class Acoustic : public Modelparameter<ValueType>
        {
        public:
            
            //! Default constructor.
            Acoustic(){};
            
            //! Destructor, releases all allocated resources.
            ~Acoustic(){};
            
            explicit Acoustic(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            explicit Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const);
            explicit Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho, IndexType partitionedIn);
            explicit Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);
            
            //! Copy Constructor.
            Acoustic(const Acoustic& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const);
            void init(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist) override;
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) override;
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho, IndexType partitionedIn);
            
            void initVelocities(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);
            
            void write(std::string filenamePWaveModulus, std::string filenamedensity, IndexType partitionedOut);
            void write(std::string filename, IndexType partitionedOut) const override;
            
            /* Getter methods for not requiered parameters */
            lama::Vector const& getSWaveModulus() override;
            lama::Vector const& getVelocityS() override;
            lama::Vector const& getTauP() override;
            lama::Vector const& getTauS() override;
            lama::Vector const& getSWaveModulusAverageXY() override;
            lama::Vector const& getSWaveModulusAverageXZ() override;
            lama::Vector const& getSWaveModulusAverageYZ() override;
            lama::Vector const& getTauSAverageXY() override;
            lama::Vector const& getTauSAverageXZ() override;
            lama::Vector const& getTauSAverageYZ() override;
            IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;
        
            void switch2velocity() override;
            void switch2modulus() override;
            
            void prepareForModelling(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, dmemo::CommunicatorPtr comm) override;
            
            /* Overloading Operators */
            KITGPI::Modelparameter::Acoustic<ValueType> operator*(lama::Scalar rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator*=(lama::Scalar rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator+(KITGPI::Modelparameter::Acoustic<ValueType> rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator+=(KITGPI::Modelparameter::Acoustic<ValueType> rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator-(KITGPI::Modelparameter::Acoustic<ValueType> rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator-=(KITGPI::Modelparameter::Acoustic<ValueType> rhs);
        
        private:
                        
            void refreshModule() override;
            void refreshVelocity() override;
            void calculateAveraging() override;
            
            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagModulus;
            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::dirtyFlagVelocity;
            using Modelparameter<ValueType>::parametrisation;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;
            
            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm ) override;
            
            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration config, dmemo::CommunicatorPtr comm );
            
            using Modelparameter<ValueType>::DensityAverageMatrixX;
            using Modelparameter<ValueType>::DensityAverageMatrixY;
            using Modelparameter<ValueType>::DensityAverageMatrixZ;
            
            using Modelparameter<ValueType>::inverseDensityAverageX;
            using Modelparameter<ValueType>::inverseDensityAverageY;
            using Modelparameter<ValueType>::inverseDensityAverageZ;
            
            /* Not requiered parameters */
            using Modelparameter<ValueType>::velocityS;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::tauP;
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::sWaveModulusAverageXY;
            using Modelparameter<ValueType>::sWaveModulusAverageXZ;
            using Modelparameter<ValueType>::sWaveModulusAverageYZ;
            using Modelparameter<ValueType>::tauSAverageXY;
            using Modelparameter<ValueType>::tauSAverageXZ;
            using Modelparameter<ValueType>::tauSAverageYZ;
            
        };
    }
}

/*! \brief Prepare modellparameter for modelling
 *
 * Refreshes the module if parameterisation is in terms of velocities
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 \param comm Communicator pointer
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::prepareForModelling(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, dmemo::CommunicatorPtr comm){
    refreshModule();
    initializeMatrices(dist,ctx,config,comm);
    this->getInverseDensity();
    calculateAveraging();
}

/*! \brief Switch the default parameterization of this class to modulus
 *
 * This will recalulcate the modulus vectors from the velocity vectors.
 * Moreover, the parametrisation value will be set to zero.
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::switch2modulus(){
    if(parametrisation==1){
        this->calcModuleFromVelocity(velocityP,density,pWaveModulus);
        dirtyFlagAveraging = true;
        dirtyFlagModulus=false;
        dirtyFlagVelocity=false;
        parametrisation=0;
    }
};

/*! \brief Switch the default parameterization of this class to velocity
 *
 * This will recalulcate the velocity vectors from the modulus vectors.
 * Moreover, the parametrisation value will be set to one.
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::switch2velocity(){
    if(parametrisation==0){
        this->calcVelocityFromModule(pWaveModulus,density,velocityP);
        dirtyFlagModulus=false;
        dirtyFlagVelocity=false;
        parametrisation=1;
    }
};

/*! \brief Refresh the velocity vectors with the module values
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::refreshVelocity(){
    if(parametrisation==0){
        this->calcVelocityFromModule(pWaveModulus,density,velocityP);
        dirtyFlagVelocity=false;
    }
};

/*! \brief Refresh the module vectors with the velocity values
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::refreshModule(){
    if(parametrisation==1){
        this->calcModuleFromVelocity(velocityP,density,pWaveModulus);
        dirtyFlagAveraging = true;
        dirtyFlagModulus=false;
    }
};


/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    init(config,ctx,dist);
}

/*! \brief Initialisation that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(Configuration::Configuration const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    if(config.get<IndexType>("ModelRead")){
        switch (config.get<IndexType>("ModelParametrisation")) {
            case 1:
                init(ctx,dist,config.get<std::string>("ModelFilename"),config.get<IndexType>("PartitionedIn"));
                break;
            case 2:
                parametrisation=1;
                initVelocities(ctx,dist,config.get<std::string>("ModelFilename"),config.get<IndexType>("PartitionedIn"));
                break;
            default:
                COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
                break;
        }
    } else {
        ValueType getPWaveModulus = config.get<ValueType>("rho") * config.get<ValueType>("velocityP")* config.get<ValueType>("velocityP");
        init(ctx,dist,getPWaveModulus,config.get<ValueType>("rho"));
    }
    
    if(config.get<IndexType>("ModelWrite")){
        write(config.get<std::string>("ModelFilename")+".out",config.get<IndexType>("PartitionedOut"));
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const)
{
    init(ctx,dist,pWaveModulus_const,rho_const);
}


/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param rho_const Density given as Scalar
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const)
{
    parametrisation=0;
    this->initModelparameter(pWaveModulus,ctx,dist,pWaveModulus_const);
    this->initModelparameter(density,ctx,dist,rho_const);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamePWaveModulus Name of file that will be read for the P-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 \param partitionedIn Partitioned input
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho, IndexType partitionedIn)
{
    init(ctx,dist,filenamePWaveModulus,filenamerho,partitionedIn);
}


/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamePWaveModulus Name of file that will be read for the P-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 \param partitionedIn Partitioned Inpiut
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho, IndexType partitionedIn)
{
    parametrisation=0;
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus,partitionedIn);
    this->initModelparameter(density,ctx,dist,filenamerho,partitionedIn);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    init(ctx,dist,filename,partitionedIn);
}


/*! \brief Initialisator that is reading Lame-models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    parametrisation=0;
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus,partitionedIn);
    this->initModelparameter(density,ctx,dist,filenamedensity,partitionedIn);
}


//! \brief Copy constructor
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType>::Acoustic(const Acoustic& rhs)
{
    pWaveModulus=rhs.pWaveModulus;
    velocityP=rhs.velocityP;
    inverseDensity=rhs.inverseDensity;
    density=rhs.density;
    dirtyFlagInverseDensity=rhs.dirtyFlagInverseDensity;
    dirtyFlagModulus=rhs.dirtyFlagModulus;
    dirtyFlagVelocity=rhs.dirtyFlagVelocity;
    parametrisation=rhs.parametrisation;
}

/*! \brief Initialisator that is reading Velocity-Vector
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Velocity-Vector "filename".vp.mtx" is added and for density "filename+".density.mtx" is added.
 \param partitionedIn Partitioned input
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::initVelocities(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    
    parametrisation=1;
    std::string filenameVelocityP=filename+".vp.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->initModelparameter(velocityP,ctx,dist,filenameVelocityP,partitionedIn);
    this->initModelparameter(density,ctx,dist,filenamedensity,partitionedIn);
    
}


/*! \brief Write model to an external file
 *
 \param filenamePWaveModulus Filename for P-wave modulus model
 \param filenamedensity Filename for Density model
 \param partitionedOut Partitioned output
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::write( std::string filenamePWaveModulus, std::string filenamedensity, IndexType partitionedOut)
{
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus,partitionedOut);
    this->writeModelparameter(density,filenamedensity,partitionedOut);
};


/*! \brief Write model to an external file
 *
 \param filename Filename to write files. For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 \param partitionedOut Partitioned output
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::write(std::string filename, IndexType partitionedOut) const
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus,partitionedOut);
    this->writeModelparameter(density,filenamedensity,partitionedOut);
};

//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration config, dmemo::CommunicatorPtr comm )
{
    initializeMatrices(dist,ctx,config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"), config.get<ValueType>("DT"), comm);
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
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, dmemo::CommunicatorPtr /*comm*/ )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    
    this->calcDensityAverageMatrixX(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixY(NX, NY, NZ, dist);
    this->calcDensityAverageMatrixZ(NX, NY, NZ, dist);
    
    DensityAverageMatrixX.setContextPtr( ctx );
    DensityAverageMatrixY.setContextPtr( ctx );
    DensityAverageMatrixZ.setContextPtr( ctx );
    
    
}

/*! \brief calculate averaged vectors
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Acoustic<ValueType>::calculateAveraging(){
    this->calculateInverseAveragedDensity(density,inverseDensityAverageX,DensityAverageMatrixX);
    this->calculateInverseAveragedDensity(density,inverseDensityAverageY,DensityAverageMatrixY);
    this->calculateInverseAveragedDensity(density,inverseDensityAverageZ,DensityAverageMatrixZ);
    dirtyFlagAveraging = false;
}


/*! \brief Get reference to S-wave modulus
 *
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulus(){
    COMMON_THROWEXCEPTION("S-wave modulus is not set for acoustic modelling")
    return(sWaveModulus);
}


/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getVelocityS(){
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return(velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getTauP(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauP);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getTauS(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauS);
}


/*! \brief Getter method for relaxation frequency */
template<typename ValueType>
ValueType KITGPI::Modelparameter::Acoustic<ValueType>::getRelaxationFrequency() const
{
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return(relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Acoustic<ValueType>::getNumRelaxationMechanisms() const
{
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return(numRelaxationMechanisms);
}


/*! \brief Get reference to S-wave modulus in xy-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXY(){
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return(sWaveModulusAverageXY);
}

/*! \brief Get reference to S-wave modulus in xz-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageXZ(){
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return(sWaveModulusAverageXZ);
}

/*! \brief Get reference to S-wave modulus in yz-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getSWaveModulusAverageYZ(){
    COMMON_THROWEXCEPTION("The averaged S-wave modulus is not set for acoustic modelling")
    return(sWaveModulusAverageYZ);
}

/*! \brief Get reference to tauS xy-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXY(){
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return(tauSAverageXY);
}

/*! \brief Get reference to tauS xz-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageXZ(){
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return(tauSAverageXZ);
}

/*! \brief Get reference to tauS yz-plane
 */
template<typename ValueType>
lama::Vector const& KITGPI::Modelparameter::Acoustic<ValueType>::getTauSAverageYZ(){
    COMMON_THROWEXCEPTION("There is no averaged tau parameter in an elastic modelling")
    return(tauSAverageYZ);
}


/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator*(lama::Scalar rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result;
    result.density = this->density * rhs;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus * rhs;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP * rhs;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}


/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> operator*(lama::Scalar lhs, KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    return rhs * lhs;
}


/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator*=(lama::Scalar rhs)
{
    return this * rhs;
}


/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator+(KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result;
    result.density = this->density + rhs.density;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus + rhs.pWaveModulus;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP + rhs.velocityP;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}


/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator+=(KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    return this + rhs;
}


/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator-(KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    KITGPI::Modelparameter::Acoustic<ValueType> result;
    result.density = this->density - rhs.density;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus - rhs.pWaveModulus;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP - rhs.velocityP;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown parametrisation! ");
    }
}


/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::Acoustic<ValueType> KITGPI::Modelparameter::Acoustic<ValueType>::operator-=(KITGPI::Modelparameter::Acoustic<ValueType> rhs)
{
    return this - rhs;
}
