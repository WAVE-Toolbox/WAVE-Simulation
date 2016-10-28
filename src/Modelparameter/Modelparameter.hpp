

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

namespace KITGPI {
    
    //! \brief Modelparameter namespace
    namespace Modelparameter {
        
        //! Abstract class for a single Modelparameter (Subsurface properties)
        /*!
         * This class handels a single modelparameter.
         * As this class is an abstract class, all constructors are protected.
         */
        template<typename ValueType>
        class Modelparameter
        {
            
        public:
            
            //! Default constructor.
            Modelparameter():dirtyFlagInverseDensity(1){};
            
            //! Default destructor.
            ~Modelparameter(){};
            
            /*! \brief Abstract initialisation function
             *
             * Standard initialisation function
             *
             \param ctx Context
             \param dist Distribution
             \param filename filename to read modelparameters (endings will be added by derived classes)
             */
            virtual void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)=0;
            
            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write modelparameters (endings will be added by derived classes)
             */
            virtual void write(std::string filename)=0;
            
            virtual lama::DenseVector<ValueType>& getDensity();
            virtual lama::DenseVector<ValueType>& getInverseDensity();
            virtual lama::DenseVector<ValueType>& getLambda();
            virtual lama::DenseVector<ValueType>& getMu();
            virtual lama::DenseVector<ValueType>& getVelocityP();
            virtual lama::DenseVector<ValueType>& getVelocityS();
            
        protected:
            
            IndexType dirtyFlagInverseDensity; //!< ==1 if inverseDensity has to be recalulated; ==0 if inverseDensity is up to date
            
            lama::DenseVector<ValueType> lambda; //!< Vector storing first Lame-Parameter.
            lama::DenseVector<ValueType> mu; //!< Vector storing first Lame-Parameter.
            lama::DenseVector<ValueType> density; //!< Vector storing Density.
            lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.
            
            lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.
            
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value);
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
            
            void calculateLame(lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity);
            void calculateLame(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity);
            
            void calculateLambda(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameDensity);
            void calculateLambda(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameVS, std::string filenameDensity);
            void calculateMu(lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecMu, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVS, std::string filenameDensity);
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
        };
    }
}

/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value)
{
    
    allocateModelparameter(vector,ctx,dist);
    
    vector.assign(value);
    
}


/*! \brief Init a single modelparameter by reading a model from an external file
 *
 *  Reads a single model from an external mtx file.
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    
    allocateModelparameter(vector,ctx,dist);
    
    readModelparameter(vector,filename);
    
    vector.redistribute(dist);
    
}


/*! \brief Write singe modelparameter to an external file
 *
 *  Write a single model to an external file.
 \param vector Single modelparameter which will be written to filename
 \param filename Name of file in which modelparameter will be written
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.writeToFile(filename);
};


/*! \brief Read a modelparameter from file
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.readFromFile(filename);
};


/*! \brief Allocate a single modelparameter
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Calculate Acoustic Lame-Vector from p-Velocity-Vector
 *  Acoustic:   lambda = rho * vP^2
 *  Elastic:    mu = rho * vS^2
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLame(lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorLame, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity)
{
    allocateModelparameter(vecV,ctx,dist);
    allocateModelparameter(vecDensity,ctx,dist);
    allocateModelparameter(vectorLame,ctx,dist);
    
    readModelparameter(vecV,filename);
    readModelparameter(vecDensity,filenameDensity);
    
    vecV.redistribute(dist);
    vecDensity.redistribute(dist);
    
    vectorLame=vecV;
    vectorLame.scale(vecV);
    vectorLame.scale(vecDensity);
    
};



/*! \brief Calculate Elastic Lame-Vector (Lambda) from VelocityP and VelocityS
 *  Elastic:    lambda = rho * (vP^2 - 2 * vS^2)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLame(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity)
{
    allocateModelparameter(vecVP,ctx,dist);
    allocateModelparameter(vecVS,ctx,dist);
    allocateModelparameter(vecDensity,ctx,dist);
    allocateModelparameter(vecLambda,ctx,dist);
    
    readModelparameter(vecVP,filename);
    readModelparameter(vecVS,filenameS);
    readModelparameter(vecDensity,filenameDensity);
    
    vecVP.redistribute(dist);
    vecVS.redistribute(dist);
    vecDensity.redistribute(dist);
    
    vecLambda=vecVS;
    vecLambda=(-2.)*vecLambda.scale(vecVS);
    vecLambda.invert();
    vecLambda.scale(vecVP);
    vecLambda.invert();
    vecLambda+=vecVP;
    vecLambda.scale(vecVP);
    vecLambda.scale(vecDensity);
    
};


/*! \brief Calculate Lambda (Acoustic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLambda(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameDensity)
{
    calculateLame(vecVP,vecDensity,vecLambda,ctx,dist,filenameVP,filenameDensity);
};


/*! \brief Calculate Lambda (Elastic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLambda(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecLambda, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameVS, std::string filenameDensity)
{
    calculateLame(vecVP,vecVS,vecDensity,vecLambda,ctx,dist,filenameVP,filenameVS,filenameDensity);
    
};


/*! \brief Calculate Mu (Elastic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateMu(lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecMu, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVS, std::string filenameDensity)
{
    calculateLame(vecVS,vecDensity,vecMu,ctx,dist,filenameVS,filenameDensity);
    
};


/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity(){
    if(dirtyFlagInverseDensity==1){
        inverseDensity.assign(density);
        inverseDensity.invert();
    }
    return(inverseDensity);
}

/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity(){
    dirtyFlagInverseDensity=1; // If density will be changed, the inverse has to be refreshed if it is accessed
    return(density);
}

/*! \brief Get reference to first Lame model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getLambda(){
    return(lambda);
}

/*! \brief Get reference to second Lame Parameter mu
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getMu(){
    return(mu);
}

/*! \brief Get reference to P-wave velocity
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP(){
    return(velocityP);
}

/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS(){
    return(velocityS);
}



