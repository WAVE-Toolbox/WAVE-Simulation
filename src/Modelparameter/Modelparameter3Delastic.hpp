
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
                        
        private:
            
            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::velocityS;
            
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
    pWaveModulus=rhs.pWaveModulus.copy();
    sWaveModulus=rhs.sWaveModulus.copy();
    density=rhs.density.copy();
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


