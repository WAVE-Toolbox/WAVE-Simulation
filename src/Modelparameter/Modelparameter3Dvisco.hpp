
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
        
        //! Class for Modelparameter for 3-D visco-elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the 3-D visco-elastic finite-difference simulation.
         */
        template<typename ValueType>
        class FD3Dvisco : public Modelparameter<ValueType>
        {
        public:
            
            //! Default constructor.
            FD3Dvisco(){};
            
            //! Destructor, releases all allocated resources.
            ~FD3Dvisco(){};
            
            FD3Dvisco(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            //! Copy Constructor.
            FD3Dvisco(const FD3Dvisco& rhs);
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const, IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            
            void calculateWaveModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void write(std::string filename);
            
        private:
            
            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
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
KITGPI::Modelparameter::FD3Dvisco<ValueType>::FD3Dvisco(Configuration::Configuration<ValueType>& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
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
KITGPI::Modelparameter::FD3Dvisco<ValueType>::FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar  tauP_const, lama::Scalar  tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
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
void KITGPI::Modelparameter::FD3Dvisco<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  pWaveModulus_const,lama::Scalar  sWaveModulus_const, lama::Scalar  rho_const, lama::Scalar tauP_const, lama::Scalar tauS_const,IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in)
{
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
KITGPI::Modelparameter::FD3Dvisco<ValueType>::FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
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
void KITGPI::Modelparameter::FD3Dvisco<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
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
KITGPI::Modelparameter::FD3Dvisco<ValueType>::FD3Dvisco(const FD3Dvisco& rhs)
{
    pWaveModulus=rhs.pWaveModulus;
    sWaveModulus=rhs.sWaveModulus;
    density=rhs.density;
    tauS=rhs.tauS;
    tauP=rhs.tauP;
    relaxationFrequency=rhs.relaxationFrequency;
    numRelaxationMechanisms=rhs.numRelaxationMechanisms;
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
void KITGPI::Modelparameter::FD3Dvisco<ValueType>::calculateWaveModulus(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenameVelocityP=filename+".vp.mtx";
    std::string filenameVelocityS=filename+".vs.mtx";
    std::string filenamedensity=filename+".density.mtx";
    std::string filenameTauP=filename+".tauP.mtx";
    std::string filenameTauS=filename+".tauS.mtx";
    
    this->calculatePWaveModulus(velocityP,density,pWaveModulus,ctx,dist,filenameVelocityP,filenamedensity);
    this->calculateSWaveModulus(velocityS,density,sWaveModulus,ctx,dist,filenameVelocityS,filenamedensity);
    this->initModelparameter(density,ctx,dist,filenamedensity);
    this->initModelparameter(tauS,ctx,dist,filenameTauS);
    this->initModelparameter(tauP,ctx,dist,filenameTauP);
}


/*! \brief Write model to an external file
 *
 \param filename For the P-wave modulus ".pWaveModulus.mtx" is added, for the second ".sWaveModulus.mtx", for density ".density.mtx", for tauP ".tauP.mtx"  and for tauS ".tauS.mtx" is added.
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dvisco<ValueType>::write(std::string filename)
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

/*! \brief Initialisation the relaxation mechanisms
 *
 \param numRelaxationMechanisms_in Number of relaxation mechanisms
 \param relaxationFrequency_in Relaxation frequency
 */
template<typename ValueType>
void KITGPI::Modelparameter::FD3Dvisco<ValueType>::initRelaxationMechanisms(IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in){
    if(numRelaxationMechanisms_in<1){
        COMMON_THROWEXCEPTION("The number of relaxation mechanisms should be >0 in an visco-elastic simulation")
    }
    if(relaxationFrequency_in<=0){
        COMMON_THROWEXCEPTION("The relaxation frequency should be >=0 in an visco-elastic simulation")
    }
    numRelaxationMechanisms=numRelaxationMechanisms_in;
    relaxationFrequency_in=relaxationFrequency;
}

