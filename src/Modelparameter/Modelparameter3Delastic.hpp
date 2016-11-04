
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

namespace KITGPI {
    
    //! \brief Modelparameter namespace
    namespace Modelparameter {
        
        //! Class for Modelparameter for 3-D elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the 3-D elastic finite-difference simulation.
         */
        template<typename ValueType>
        class FD3Delastic : public Modelparameter<ValueType>
        {
        public:
            
            //! Default constructor.
            FD3Delastic(){};
            
            //! Destructor, releases all allocated resources.
            ~FD3Delastic(){};
            
            FD3Delastic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus,std::string filenameSWaveModulus, std::string filenamerho);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            //! Copy Constructor.
            FD3Delastic(const FD3Delastic& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus,lama::Scalar  sWaveModulus, lama::Scalar  rho);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus,std::string filenameSWaveModulus, std::string filenamerho);
            
            void calculateWaveModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void write(std::string filenamePWaveModulus,std::string filenameSWaveModulus, std::string filenamedensity);
            void write(std::string filename);
            
            /* Getter methods for not requiered parameters */
            lama::DenseVector<ValueType>& getTauP();
            lama::DenseVector<ValueType>& getTauS();
            IndexType getNumRelaxationMechanisms();
            ValueType getRelaxationFrequency();
            
            /* Overloading Operators */
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator*(lama::Scalar rhs);
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator*=(lama::Scalar rhs);
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator+(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator+=(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator-(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Delastic<ValueType> operator-=(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs);

        private:
            
            /*! \brief Prepare the model parameters for modelling */
            /* Nothing has to be done here */
            void prepareForModelling(){};
            
            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagParametrisation;
            using Modelparameter<ValueType>::Parametrisation;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::velocityS;
            
            /* Not requiered parameters */
            using Modelparameter<ValueType>::tauP;
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            
        };
    }
}

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    if(config.getModelRead()){
        switch (config.getModelParametrisation()) {
            case 1:
                init(ctx,dist,config.getModelFilename());
                break;
            case 2:
                dirtyFlagParametrisation=1;
                calculateWaveModulus(ctx,dist,config.getModelFilename());
                break;
            default:
                COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
                break;
        }
    } else {
        init(ctx,dist,config.getPWaveModulus(),config.getSWaveModulus(),config.getRho());
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
 \param rho Density given as Scalar
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho)
{
    init(ctx,dist,pWaveModulus_const,sWaveModulus_const,rho);
}


/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param pWaveModulus_const P-wave modulus given as Scalar
 \param sWaveModulus_const S-wave modulus given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho)
{
    this->initModelparameter(pWaveModulus,ctx,dist,pWaveModulus_const);
    this->initModelparameter(sWaveModulus,ctx,dist,sWaveModulus_const);
    this->initModelparameter(density,ctx,dist,rho);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamePWaveModulus Name of file that will be read for the P-wave modulus.
 \param filenameSWaveModulus Name of file that will be read for the S-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus,std::string filenameSWaveModulus, std::string filenamerho)
{
    init(ctx,dist,filenamePWaveModulus,filenameSWaveModulus,filenamerho);
}


/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamePWaveModulus Name of file that will be read for the P-wave modulus.
 \param filenameSWaveModulus Name of file that will be read for the S-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenameSWaveModulus, std::string filenamerho)
{
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus);
    this->initModelparameter(density,ctx,dist,filenamerho);
    this->initModelparameter(sWaveModulus,ctx,dist,filenameSWaveModulus);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    init(ctx,dist,filename);
}


/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenameSWaveModulus=filename+".sWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus);
    this->initModelparameter(sWaveModulus,ctx,dist,filenameSWaveModulus);
    this->initModelparameter(density,ctx,dist,filenamedensity);
}


//! \brief Copy constructor
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(const FD3Delastic& rhs)
{
    pWaveModulus=rhs.pWaveModulus;
    sWaveModulus=rhs.sWaveModulus;
    velocityP=rhs.velocityP;
    velocityS=rhs.velocityS;
    density=rhs.density;
    dirtyFlagInverseDensity=rhs.dirtyFlagInverseDensity;
    dirtyFlagParametrisation=rhs.dirtyFlagParametrisation;
    Parametrisation=rhs.Parametrisation;
    inverseDensity=rhs.inverseDensity;
}

/*! \brief Initialisator that is reading Velocity-Vector from an external files and calculates pWaveModulus
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Velocity-Vector "filename".vp.mtx" is added and for density "filename+".density.mtx" is added.
 *
 *  Calculates pWaveModulus with
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::calculateWaveModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenameVelocityP=filename+".vp.mtx";
    std::string filenameVelocityS=filename+".vs.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->calculatePWaveModulus(velocityP,density,pWaveModulus,ctx,dist,filenameVelocityP,filenamedensity);
    this->calculateSWaveModulus(velocityS,density,sWaveModulus,ctx,dist,filenameVelocityS,filenamedensity);
    this->initModelparameter(density,ctx,dist,filenamedensity);
    
}


/*! \brief Write model to an external file
 *
 \param filenamePWaveModulus Filename for P-wave modulus model
 \param filenameSWaveModulus Name of file that will be read for the S-wave modulus.
 \param filenamedensity Filename for Density model
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::write( std::string filenamePWaveModulus, std::string filenameSWaveModulus, std::string filenamedensity)
{
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus);
    this->writeModelparameter(density,filenamedensity);
    this->writeModelparameter(sWaveModulus,filenameSWaveModulus);
};


/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx" and for density ".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::write(std::string filename)
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenameSWaveModulus=filename+".sWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus);
    this->writeModelparameter(sWaveModulus,filenameSWaveModulus);
    this->writeModelparameter(density,filenamedensity);
};

/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getTauP(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauP);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getTauS(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauS);
}

/*! \brief Getter method for relaxation frequency */
template<typename ValueType>
ValueType KITGPI::Modelparameter::FD3Delastic<ValueType>::getRelaxationFrequency(){
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return(relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template<typename ValueType>
IndexType KITGPI::Modelparameter::FD3Delastic<ValueType>::getNumRelaxationMechanisms(){
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return(numRelaxationMechanisms);
}


/*! \brief Overloading * Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator*(lama::Scalar rhs)
{
    KITGPI::Modelparameter::FD3Delastic<ValueType> result;
    result.density = this->density * rhs;
    if (Parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus * rhs;
        result.sWaveModulus= this->sWaveModulus * rhs;
        return result;
    } if (Parametrisation==1) {
        result.velocityP= this->velocityP * rhs;
        result.velocityS= this->velocityS * rhs;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}


/*! \brief free function to multiply
 *
 \param lhs Scalar factor with which the vectors are multiplied.
 \param rhs Vector
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> operator*(lama::Scalar lhs, KITGPI::Modelparameter::FD3Delastic<ValueType> rhs)
{
    return rhs * lhs;
}


/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator*=(lama::Scalar rhs)
{
    return this * rhs;
}


/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator+(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs)
{
    KITGPI::Modelparameter::FD3Delastic<ValueType> result;
    result.density = this->density + rhs.density;
    if (Parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus + rhs.pWaveModulus;
        result.sWaveModulus= this->sWaveModulus + rhs.sWaveModulus;
        return result;
    } if (Parametrisation==1) {
        result.velocityP= this->velocityP + rhs.velocityP;
        result.velocityS= this->velocityS + rhs.velocityS;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}


/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator+=(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs)
{
    return this + rhs;
}


/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator-(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs)
{
    KITGPI::Modelparameter::FD3Delastic<ValueType> result;
    result.density = this->density - rhs.density;
    if (Parametrisation==0) {
        result.pWaveModulus= this->pWaveModulus - rhs.pWaveModulus;
        result.sWaveModulus= this->sWaveModulus - rhs.sWaveModulus;
        return result;
    } if (Parametrisation==1) {
        result.velocityP= this->velocityP - rhs.velocityP;
        result.velocityS= this->velocityS - rhs.velocityS;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}


/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType> KITGPI::Modelparameter::FD3Delastic<ValueType>::operator-=(KITGPI::Modelparameter::FD3Delastic<ValueType> rhs)
{
    return this - rhs; 
}


