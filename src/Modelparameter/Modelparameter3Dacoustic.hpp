
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
        
        //! Class for Modelparameter for 3-D acoustic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the 3-D acoustic finite-difference simulation.
         */
        template<typename ValueType>
        class FD3Dacoustic : public Modelparameter<ValueType>
        {
        public:
            
            //! Default constructor.
            FD3Dacoustic(){};
            
            //! Destructor, releases all allocated resources.
            ~FD3Dacoustic(){};
            
            FD3Dacoustic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const);
            FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho);
            FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            //! Copy Constructor.
            FD3Dacoustic(const FD3Dacoustic& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho);
            
            void calculateModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void write(std::string filenamePWaveModulus, std::string filenamedensity);
            void write(std::string filename);
            
            /* Getter methods for not requiered parameters */
            lama::DenseVector<ValueType>& getSWaveModulus();
            lama::DenseVector<ValueType>& getVelocityS();
            lama::DenseVector<ValueType>& getTauP();
            lama::DenseVector<ValueType>& getTauS();
            IndexType getNumRelaxationMechanisms();
            ValueType getRelaxationFrequency();
            
            /* Overloading Operators */
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator*=(ValueType rhs);
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator+(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator+=(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator-(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs);
            KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator-=(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs);
            
        private:
            
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
KITGPI::Modelparameter::FD3Dacoustic<ValueType>::FD3Dacoustic(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    if(config.getModelRead()){
        switch (config.getModelParametrisation()) {
            case 1:
                init(ctx,dist,config.getModelFilename());
                break;
            case 2:
                Parametrisation=1;
                calculateModulus(ctx,dist,config.getModelFilename());
                break;
            default:
                COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
                break;
        }
    } else {
        init(ctx,dist,config.getPWaveModulus(),config.getRho());
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
 \param rho_const Density given as Scalar
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType>::FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const)
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
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const, lama::Scalar  rho_const)
{
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
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType>::FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho)
{
    init(ctx,dist,filenamePWaveModulus,filenamerho);
}


/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamePWaveModulus Name of file that will be read for the P-wave modulus.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamePWaveModulus, std::string filenamerho)
{
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus);
    this->initModelparameter(density,ctx,dist,filenamerho);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType>::FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    init(ctx,dist,filename);
}


/*! \brief Initialisator that is reading Lame-models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added and for density "filename+".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->initModelparameter(pWaveModulus,ctx,dist,filenamePWaveModulus);
    this->initModelparameter(density,ctx,dist,filenamedensity);
}


//! \brief Copy constructor
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType>::FD3Dacoustic(const FD3Dacoustic& rhs)
{
    pWaveModulus=rhs.pWaveModulus;
    density=rhs.density;
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
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::calculateModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenameVelocityP=filename+".vp.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->calculatePWaveModulus(velocityP,density,pWaveModulus,ctx,dist,filenameVelocityP,filenamedensity);
    this->initModelparameter(density,ctx,dist,filenamedensity);
    
}


/*! \brief Write model to an external file
 *
 \param filenamePWaveModulus Filename for P-wave modulus model
 \param filenamedensity Filename for Density model
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::write( std::string filenamePWaveModulus, std::string filenamedensity)
{
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus);
    this->writeModelparameter(density,filenamedensity);
};


/*! \brief Write model to an external file
 *
 \param filename Filename to write files. For the P-wave modulus ".pWaveModulus.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dacoustic<ValueType>::write(std::string filename)
{
    std::string filenamePWaveModulus=filename+".pWaveModulus.mtx";
    std::string filenamedensity=filename+".density.mtx";
    this->writeModelparameter(pWaveModulus,filenamePWaveModulus);
    this->writeModelparameter(density,filenamedensity);
};



/*! \brief Get reference to S-wave modulus
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getSWaveModulus(){
    COMMON_THROWEXCEPTION("S-wave modulus is not set for acoustic modelling")
    return(sWaveModulus);
}


/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getVelocityS(){
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an acoustic simulation.")
    return(velocityS);
}

/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getTauP(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauP);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getTauS(){
    COMMON_THROWEXCEPTION("There is no tau parameter in an elastic modelling")
    return(tauS);
}


/*! \brief Getter method for relaxation frequency */
template<typename ValueType>
ValueType KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getRelaxationFrequency(){
    COMMON_THROWEXCEPTION("There is no relaxationFrequency parameter in an elastic modelling")
    return(relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template<typename ValueType>
IndexType KITGPI::Modelparameter::FD3Dacoustic<ValueType>::getNumRelaxationMechanisms(){
    COMMON_THROWEXCEPTION("There is no numRelaxationMechanisms parameter in an elastic modelling")
    return(numRelaxationMechanisms);
}


/*! \brief Overloading * Operation 
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator*(ValueType rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus * rhs;
        result.density = this->density * rhs;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP * rhs;
        result.density = this->density * rhs;
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
KITGPI::Modelparameter::FD3Dacoustic<ValueType> operator*(ValueType lhs, KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs)
{
    return rhs * lhs;
}


/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator*=(ValueType rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus * rhs;
        result.density = this->density * rhs;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP * rhs;
        result.density = this->density * rhs;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}


/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator+(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus + rhs.pWaveModulus;
        result.density = this->density + rhs.density;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP + rhs.velocityP;
        result.density = this->density + rhs.density;
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
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator+=(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus + rhs.pWaveModulus;
        result.density = this->density + rhs.density;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP + rhs.velocityP;
        result.density = this->density + rhs.density;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}


/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator-(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus - rhs.pWaveModulus;
        result.density = this->density - rhs.density;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP - rhs.velocityP;
        result.density = this->density - rhs.density;
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
KITGPI::Modelparameter::FD3Dacoustic<ValueType> KITGPI::Modelparameter::FD3Dacoustic<ValueType>::operator-=(KITGPI::Modelparameter::FD3Dacoustic<ValueType> rhs)
{
    if (Parametrisation==0) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.pWaveModulus= this->pWaveModulus - rhs.pWaveModulus;
        result.density = this->density - rhs.density;
        return result;
    } if (Parametrisation==1) {
        KITGPI::Modelparameter::FD3Dacoustic<ValueType> result;
        result.velocityP= this->velocityP - rhs.velocityP;
        result.density = this->density - rhs.density;
        return result;
    } else {
        COMMON_THROWEXCEPTION(" Unknown Parametrisation! ");
    }
}
