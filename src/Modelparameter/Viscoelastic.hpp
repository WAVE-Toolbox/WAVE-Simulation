
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include "../Common/HostPrint.hpp"

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "Modelparameter.hpp"

namespace KITGPI {
    
    //! \brief Modelparameter namespace
    namespace Modelparameter {
        
        //! Class for Modelparameter for visco-elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the visco-elastic finite-difference simulation.
         */
        template<typename ValueType>
        class Viscoelastic : public Modelparameter<ValueType>
        {
        public:
            
            //! Default constructor.
            Viscoelastic(){};
            
            //! Destructor, releases all allocated resources.
            ~Viscoelastic(){};
            
            Viscoelastic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            Viscoelastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            Viscoelastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            //! Copy Constructor.
            Viscoelastic(const Viscoelastic& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            
            void initVelocities(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void write(std::string filename);
            
            void prepareForModelling();
            
            void switch2velocity();
            void switch2modulus();
            
            /* Overloading Operators */
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator*(lama::Scalar rhs);
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator*=(lama::Scalar rhs);
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator+(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs);
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator+=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs);
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator-(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs);
            KITGPI::Modelparameter::Viscoelastic<ValueType> operator-=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs);
            
        private:
            
            void refreshModule();
            void refreshVelocity();
            
            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagModulus;
            using Modelparameter<ValueType>::dirtyFlagVelocity;
            using Modelparameter<ValueType>::parametrisation;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::velocityS;
            
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::tauP;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            
            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm );
            void runAveraging(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx);
            
            using Modelparameter<ValueType>::xAvDensityMatrix;
            using Modelparameter<ValueType>::yAvDensityMatrix;
            using Modelparameter<ValueType>::zAvDensityMatrix;
            using Modelparameter<ValueType>::xyAvSWaveModulusMatrix;
            using Modelparameter<ValueType>::xzAvSWaveModulusMatrix;
            using Modelparameter<ValueType>::yzAvSWaveModulusMatrix;
            
            using Modelparameter<ValueType>::xAvDensity;
            using Modelparameter<ValueType>::yAvDensity;
            using Modelparameter<ValueType>::zAvDensity;
            using Modelparameter<ValueType>::xyAvSWaveModulus;
            using Modelparameter<ValueType>::xzAvSWaveModulus;
            using Modelparameter<ValueType>::yzAvSWaveModulus;
            using Modelparameter<ValueType>::xyAvTauS;
            using Modelparameter<ValueType>::xzAvTauS;
            using Modelparameter<ValueType>::yzAvTauS;
            
            
        };
    }
}


/*! \brief Switch the default parameterization of this class to modulus
 *
 * This will recalulcate the modulus vectors from the velocity vectors.
 * Moreover, the parametrisation value will be set to zero.
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::switch2modulus(){
    if(parametrisation==1){
        this->calcModuleFromVelocity(velocityP,density,pWaveModulus);
        this->calcModuleFromVelocity(velocityS,density,sWaveModulus);
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
void KITGPI::Modelparameter::Viscoelastic<ValueType>::switch2velocity(){
    if(parametrisation==0){
        this->calcVelocityFromModule(pWaveModulus,density,velocityP);
        this->calcVelocityFromModule(sWaveModulus,density,velocityS);
        dirtyFlagModulus=false;
        dirtyFlagVelocity=false;
        parametrisation=1;
    }
};

/*! \brief Refresh the velocity vectors with the module values
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::refreshVelocity(){
    if(parametrisation==0){
        this->calcVelocityFromModule(pWaveModulus,density,velocityP);
        this->calcVelocityFromModule(sWaveModulus,density,velocityS);
        dirtyFlagVelocity=false;
    }
};

/*! \brief Refresh the module vectors with the velocity values
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::refreshModule(){
    if(parametrisation==1){
        this->calcModuleFromVelocity(velocityP,density,pWaveModulus);
        this->calcModuleFromVelocity(velocityS,density,sWaveModulus);
        dirtyFlagModulus=false;
    }
};


/*! \brief Prepare modellparameter for visco-elastic modelling
 *
 * Applies Equation 12 from Bohlen 2002 and refreshes the module
 *
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::prepareForModelling(){
    
    refreshModule();
    
    /* Set circular frequency w = 2 * pi * relaxation frequency */
    ValueType w_ref = 2.0 * M_PI * relaxationFrequency;
    ValueType tauSigma= 1.0 / ( 2.0 * M_PI * relaxationFrequency );
    
    ValueType sum= w_ref*w_ref*tauSigma*tauSigma / (1.0 + w_ref*w_ref*tauSigma*tauSigma );
    
    lama::DenseVector<ValueType> temp(tauS.getDistributionPtr());
    
    /* Scaling the S-wave Modulus */
    temp = 1.0;
    temp += sum * tauS;
    temp.invert();
    sWaveModulus.scale(temp);
    
    /* Scaling the P-wave Modulus */
    temp = 1.0;
    temp += sum * tauP;
    temp.invert();
    pWaveModulus.scale(temp);
    
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    if(config.getModelRead()){
        switch (config.getModelParametrisation()) {
            case 1:
                init(ctx,dist,config.getModelFilename());
                break;
            case 2:
                initVelocities(ctx,dist,config.getModelFilename());
                break;
            default:
                COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
                break;
        }
        initRelaxationMechanisms(config.getNumRelaxationMechanisms(), config.getRelaxationFrequency());
        
    } else {
        init(ctx,dist,config.getPWaveModulus(),config.getSWaveModulus(),config.getRho(),config.getTauP(),config.getTauS(),config.getNumRelaxationMechanisms(), config.getRelaxationFrequency());
    }
    
    if(config.getModelWrite()){
        write(config.getModelFilename()+".out");
    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
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
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar  tauP_const, lama::Scalar  tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    init(ctx,dist,pWaveModulus_const,sWaveModulus_const,rho_const,tauP_const,tauS_const);
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
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
    parametrisation=0;
    initRelaxationMechanisms(numRelaxationMechanisms_in, relaxationFrequency_in);
    this->initModelparameter(pWaveModulus,ctx,dist,pWaveModulus_const);
    this->initModelparameter(sWaveModulus,ctx,dist,sWaveModulus_const);
    this->initModelparameter(density,ctx,dist,rho_const);
    this->initModelparameter(tauS,ctx,dist,tauS_const);
    this->initModelparameter(tauP,ctx,dist,tauP_const);
    
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    init(ctx,dist,filename);
}


/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    parametrisation=0;
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenameSWaveModulus=filename+".sWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    std::string filenameTauP=filename+".tauP.mtx";
    std::string filenameTauS=filename+".tauS.mtx";
    
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus);
    this->initModelparameter(sWaveModulus,ctx,dist,filenameSWaveModulus);
    this->initModelparameter(density,ctx,dist,filenamedensity);
    this->initModelparameter(tauS,ctx,dist,filenameTauS);
    this->initModelparameter(tauP,ctx,dist,filenameTauP);
    
}


//! \brief Copy constructor
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType>::Viscoelastic(const Viscoelastic& rhs)
{
    pWaveModulus=rhs.pWaveModulus;
    sWaveModulus=rhs.sWaveModulus;
    velocityP=rhs.velocityP;
    velocityS=rhs.velocityS;
    density=rhs.density;
    tauS=rhs.tauS;
    tauP=rhs.tauP;
    relaxationFrequency=rhs.relaxationFrequency;
    numRelaxationMechanisms=rhs.numRelaxationMechanisms;
    dirtyFlagInverseDensity=rhs.dirtyFlagInverseDensity;
    dirtyFlagModulus=rhs.dirtyFlagModulus;
    dirtyFlagVelocity=rhs.dirtyFlagVelocity;
    parametrisation=rhs.parametrisation;
    inverseDensity=rhs.inverseDensity;
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 *
 *  Calculates pWaveModulus with
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initVelocities(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    parametrisation=1;
    std::string filenameVelocityP=filename+".vp.mtx";
    std::string filenameVelocityS=filename+".vs.mtx";
    std::string filenamedensity=filename+".density.mtx";
    std::string filenameTauP=filename+".tauP.mtx";
    std::string filenameTauS=filename+".tauS.mtx";
    
    this->initModelparameter(velocityS,ctx,dist,filenameVelocityS);
    this->initModelparameter(velocityP,ctx,dist,filenameVelocityP);
    this->initModelparameter(density,ctx,dist,filenamedensity);
    this->initModelparameter(tauS,ctx,dist,filenameTauS);
    this->initModelparameter(tauP,ctx,dist,filenameTauP);
}


/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::write(std::string filename)
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenameSWaveModulus=filename+".sWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    std::string filenameTauP=filename+".tauP.mtx";
    std::string filenameTauS=filename+".tauS.mtx";
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus);
    this->writeModelparameter(sWaveModulus,filenameSWaveModulus);
    this->writeModelparameter(density,filenamedensity);
    this->writeModelparameter(tauP,filenameTauP);
    this->writeModelparameter(tauS,filenameTauS);
    
};

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
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType /*DH*/, ValueType /*DT*/, dmemo::CommunicatorPtr comm )
{
    
    SCAI_REGION( "initializeMatrices" )
    
    HOST_PRINT( comm, "Initialization of the averaging matrices." );
    
    this->calc_xAvDensityMatrix(NX, NY, NZ, dist);
    this->calc_yAvDensityMatrix(NX, NY, NZ, dist);
    this->calc_zAvDensityMatrix(NX, NY, NZ, dist);
    this->calc_xyAvSWaveModulusMatrix(NX, NY, NZ, dist);
    this->calc_xzAvSWaveModulusMatrix(NX, NY, NZ, dist);
    this->calc_yzAvSWaveModulusMatrix(NX, NY, NZ, dist);
    
    xAvDensityMatrix.setContextPtr( ctx );
    yAvDensityMatrix.setContextPtr( ctx );
    zAvDensityMatrix.setContextPtr( ctx );
    xyAvSWaveModulusMatrix.setContextPtr( ctx );
    xzAvSWaveModulusMatrix.setContextPtr( ctx );
    yzAvSWaveModulusMatrix.setContextPtr( ctx );
    
    HOST_PRINT( comm, "Finished with initialization of the matrices!\n" );
    
    
}


/*! \brief calculate averaged vectors
 *
 \param dist Distribution
 \param ctx Context
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::runAveraging(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
    this->calculateAveragedDensity(xAvDensity,xAvDensityMatrix,ctx,dist);
    this->calculateAveragedDensity(yAvDensity,yAvDensityMatrix,ctx,dist);
    this->calculateAveragedDensity(zAvDensity,zAvDensityMatrix,ctx,dist);
    
    this->calculateAveragedSWaveModulus(xyAvSWaveModulus,xyAvSWaveModulusMatrix,ctx,dist);
    this->calculateAveragedSWaveModulus(xzAvSWaveModulus,xzAvSWaveModulusMatrix,ctx,dist);
    this->calculateAveragedSWaveModulus(yzAvSWaveModulus,yzAvSWaveModulusMatrix,ctx,dist);
    
    this->calculateAveragedTauS(xyAvTauS,xyAvSWaveModulusMatrix,ctx,dist);
    this->calculateAveragedTauS(xzAvTauS,xzAvSWaveModulusMatrix,ctx,dist);
    this->calculateAveragedTauS(yzAvTauS,yzAvSWaveModulusMatrix,ctx,dist);
}


/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template<typename ValueType>
void KITGPI::Modelparameter::Viscoelastic<ValueType>::initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in){
    if(numRelaxationMechanisms_in<1){
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-elastic simulation")
    }
    if(relaxationFrequency_in<=0){
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-elastic simulation")
    }
    numRelaxationMechanisms=numRelaxationMechanisms_in;
    relaxationFrequency_in=relaxationFrequency;
}


/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*(lama::Scalar rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density * rhs;
    result.tauS = this->tauS * rhs;
    result.tauP = this->tauP * rhs;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus * rhs;
        result.sWaveModulus= this->sWaveModulus * rhs;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP * rhs;
        result.velocityS= this->velocityS * rhs;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> operator*(lama::Scalar lhs, KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return rhs * lhs;
}


/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator*=(lama::Scalar rhs)
{
    return this * rhs;
}


/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density + rhs.density;
    result.tauS = this->tauS + rhs.tauS;
    result.tauP = this->tauP + rhs.tauP;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus + rhs.pWaveModulus;
        result.sWaveModulus= this->sWaveModulus + rhs.sWaveModulus;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP + rhs.velocityP;
        result.velocityS= this->velocityS + rhs.velocityS;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator+=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return this + rhs;
}


/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    KITGPI::Modelparameter::Viscoelastic<ValueType> result;
    result.density = this->density - rhs.density;
    result.tauS = this->tauS - rhs.tauS;
    result.tauP = this->tauP - rhs.tauP;
    if (parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus - rhs.pWaveModulus;
        result.sWaveModulus= this->sWaveModulus - rhs.sWaveModulus;
        return result;
    } if (parametrisation==1) {
        result.velocityP= this->velocityP - rhs.velocityP;
        result.velocityS= this->velocityS - rhs.velocityS;
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
KITGPI::Modelparameter::Viscoelastic<ValueType> KITGPI::Modelparameter::Viscoelastic<ValueType>::operator-=(KITGPI::Modelparameter::Viscoelastic<ValueType> rhs)
{
    return this - rhs;
}
