

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
            Modelparameter():dirtyFlagInverseDensity(1),numRelaxationMechanisms(0){};
            
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
            virtual lama::DenseVector<ValueType>& getPWaveModulus();
            virtual lama::DenseVector<ValueType>& getSWaveModulus();
            virtual lama::DenseVector<ValueType>& getVelocityP();
            virtual lama::DenseVector<ValueType>& getVelocityS();
            
            virtual lama::DenseVector<ValueType>& getTauP();
            virtual lama::DenseVector<ValueType>& getTauS();
            
            virtual IndexType getNumRelaxationMechanisms();
            virtual ValueType getRelaxationFrequency();
            
        protected:
            
            IndexType dirtyFlagInverseDensity; //!< ==1 if inverseDensity has to be recalulated; ==0 if inverseDensity is up to date
            
            lama::DenseVector<ValueType> pWaveModulus; //!< Vector storing P-wave modulus.
            lama::DenseVector<ValueType> sWaveModulus; //!< Vector storing P-wave modulus.
            lama::DenseVector<ValueType> density; //!< Vector storing Density.
            lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.
            
            lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.
            
            lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.
            
            IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency; //!< Relaxation Frequency
            
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value);
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
            
            void calculateModulus(lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity);
            
            void calculatePWaveModulus(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecPWaveModulus, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameDensity);
            void calculateSWaveModulus(lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecSWaveModulus, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVS, std::string filenameDensity);
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
        };
    }
}


/*! \brief Getter method for relaxation frequency */
template<typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getRelaxationFrequency(){
    return(relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getNumRelaxationMechanisms(){
    return(numRelaxationMechanisms);
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

/*! \brief Calculate S and P wave modulus from velocities
 *  Acoustic:   pWaveModulus = rho * vP^2
 *  Elastic:    sWaveModulus = rho * vS^2
 \param vecV Velocity-Vector which will be used in the calculation (vP: Acoustic, vS: Elastic)
 \param vecDensity Density-Vector which will be used in the calculation
 \param vectorModulus Lame-Vector which is calculated
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param filenameDensity Location of external density-file which will be read in
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateModulus(lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorModulus, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity)
{
    allocateModelparameter(vecV,ctx,dist);
    allocateModelparameter(vecDensity,ctx,dist);
    allocateModelparameter(vectorModulus,ctx,dist);
    
    readModelparameter(vecV,filename);
    readModelparameter(vecDensity,filenameDensity);
    
    vecV.redistribute(dist);
    vecDensity.redistribute(dist);
    
    vectorModulus=vecV;
    vectorModulus.scale(vecV);
    vectorModulus.scale(vecDensity);
    
};


/*! \brief Calculate P-Wave Modulus
 *
 \param vecVP Velocity-Vector (VP) which will be used to calculete pWaveModulus
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecPWaveModulus Lame-Vector which is calculated
 \param ctx Context
 \param dist Distribution
 \param filenameVP Location of external VP-file which will be read in
 \param filenameDensity Location of external density-file which will be read in
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculatePWaveModulus(lama::DenseVector<ValueType>& vecVP, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecPWaveModulus, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVP, std::string filenameDensity)
{
    calculateModulus(vecVP,vecDensity,vecPWaveModulus,ctx,dist,filenameVP,filenameDensity);
};



/*! \brief Calculate S-Wave Modulus
 *
 \param vecVS Velocity-Vector (VS) which will be used to calculete sWaveModulus
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecSWaveModulus sWaveModulus-Vector which is calculated
 \param ctx Context
 \param dist Distribution
 \param filenameVS Location of external VS-file which will be read in
 \param filenameDensity Location of external density-file which will be read in
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateSWaveModulus(lama::DenseVector<ValueType>& vecVS, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecSWaveModulus, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenameVS, std::string filenameDensity)
{
    calculateModulus(vecVS,vecDensity,vecSWaveModulus,ctx,dist,filenameVS,filenameDensity);
    
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
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus(){
    return(pWaveModulus);
}

/*! \brief Get reference to second Lame Parameter sWaveModulus
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus(){
    return(sWaveModulus);
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

/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauP(){
    return(tauP);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauS(){
    return(tauS);
}

