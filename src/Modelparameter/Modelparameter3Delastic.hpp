
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
            FD3Delastic():dirtyFlagInverseDensity(1){};
            
            //! Destructor, releases all allocated resources.
            ~FD3Delastic(){};
            
            FD3Delastic(Configuration::Configuration<ValueType> config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M_const,lama::Scalar  Mu_const, lama::Scalar  rho);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameM,std::string filenameMu, std::string filenamerho);
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            //! Copy Constructor.
            FD3Delastic(const FD3Delastic& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  lambda,lama::Scalar  mu, lama::Scalar  rho);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameM,std::string filenameMu, std::string filenamerho);
            
            void write(std::string filenameM,std::string filenameMu, std::string filenamedensity);
            void write(std::string filename);
            
            /* Getter routines for modelparameters */
            lama::DenseVector<ValueType>& getDensity();
            lama::DenseVector<ValueType>& getInverseDensity();
            lama::DenseVector<ValueType>& getLambda();
            lama::DenseVector<ValueType>& getMu();
            lama::DenseVector<ValueType>& getVelocityP();
            lama::DenseVector<ValueType>& getVelocityS();
            
        private:
            
            IndexType dirtyFlagInverseDensity; //!< ==1 if inverseDensity has to be recalulated; ==0 if inverseDensity is up to date
            
            lama::DenseVector<ValueType> lambda; //!< Vector storing first Lame-Parameter.
            lama::DenseVector<ValueType> mu; //!< Vector storing first Lame-Parameter.
            lama::DenseVector<ValueType> density; //!< Vector storing Density.
            lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.
            
            lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.
            
            
            
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
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(Configuration::Configuration<ValueType> /*config*/, hmemo::ContextPtr /*ctx*/, dmemo::DistributionPtr /*dist*/)
:dirtyFlagInverseDensity(1)
{
//    if(config.getModelRead()){
//        init(ctx,dist,config.getModelFilename());
//    } else {
//        init(ctx,dist,config.getLambda(),config.getRho());
//    }
//    
//    if(config.getModelWrite()){
//        write(config.getModelFilename()+".out");
//    }
}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param M_const First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M_const,lama::Scalar  Mu_const, lama::Scalar  rho)
:dirtyFlagInverseDensity(1)
{
    init(ctx,dist,M_const,Mu_const,rho);
}


/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param M_const First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M_const,lama::Scalar  Mu_const, lama::Scalar  rho)
{
    this->initModelparameter(lambda,ctx,dist,M_const);
    this->initModelparameter(mu,ctx,dist,Mu_const);
    this->initModelparameter(density,ctx,dist,rho);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenameM Name of file that will be read for the first Lame-parameter.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameM,std::string filenameMu, std::string filenamerho)
:dirtyFlagInverseDensity(1)
{
    init(ctx,dist,filenameM,filenameMu,filenamerho);
}


/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenameM Name of file that will be read for the first Lame-parameter.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameM, std::string filenameMu, std::string filenamerho)
{
    this->initModelparameter(lambda,ctx,dist,filenameM);
    this->initModelparameter(density,ctx,dist,filenamerho);
    this->initModelparameter(mu,ctx,dist,filenameMu);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Lame-parameter ".M.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
:dirtyFlagInverseDensity(1)
{
    init(ctx,dist,filename);
}


/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Lame-parameter ".M.mtx" is added and for density "filename+".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenameM=filename+".Lambda.mtx";
    std::string filenameMu=filename+".mu.mtx";
    std::string filenamedensity=filename+".density.mtx";
    
    this->initModelparameter(lambda,ctx,dist,filenameM);
    this->initModelparameter(mu,ctx,dist,filenameMu);
    this->initModelparameter(density,ctx,dist,filenamedensity);
}


//! \brief Copy constructor
template<typename ValueType>
KITGPI::Modelparameter::FD3Delastic<ValueType>::FD3Delastic(const FD3Delastic& rhs)
:dirtyFlagInverseDensity(1)
{
    lambda=rhs.M.copy();
    mu=rhs.mu.copy();
    density=rhs.density.copy();
}


/*! \brief Write model to an external file
 *
 \param filenameM Filename for first Lame-Parameter model
 \param filenamedensity Filename for Density model
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::write( std::string filenameM, std::string filenameMu, std::string filenamedensity)
{
    this->writeModelparameter(lambda,filenameM);
    this->writeModelparameter(density,filenamedensity);
    this->writeModelparameter(mu,filenameMu);
};


/*! \brief Write model to an external file
 *
 \param filename Filename to write files. For the first Lame-parameter ".M.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Delastic<ValueType>::write(std::string filename)
{
    std::string filenameM=filename+".Lambda.mtx";
    std::string filenameMu=filename+".mu.mtx";
    std::string filenamedensity=filename+".density.mtx";
    this->writeModelparameter(lambda,filenameM);
    this->writeModelparameter(mu,filenameM);
    this->writeModelparameter(density,filenamedensity);
};



/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getInverseDensity(){
    if(dirtyFlagInverseDensity==1){
        inverseDensity.assign(density);
        inverseDensity.invert();
    }
    return(inverseDensity);
}

/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getDensity(){
    dirtyFlagInverseDensity=1; // If density will be changed, the inverse has to be refreshed if it is accessed
    return(density);
}

/*! \brief Get reference to first Lame model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getLambda(){
    return(lambda);
}

/*! \brief Get reference to second Lame Parameter mu
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getMu(){
    return(mu);
}

/*! \brief Get reference to P-wave velocity
 *
 * Not yet implemented
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getVelocityP(){
    COMMON_THROWEXCEPTION("Conversation to P-wave velocity is not yet implemented")
    return(velocityP);
}

/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::FD3Delastic<ValueType>::getVelocityS(){
    COMMON_THROWEXCEPTION("The S-wave velocity is not defined in an elastic simulation.")
    return(velocityS);
}

