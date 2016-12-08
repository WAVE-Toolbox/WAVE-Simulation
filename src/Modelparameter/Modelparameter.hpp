

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

#include "../Partitioning/PartitionedInOut.hpp"

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
            Modelparameter():dirtyFlagInverseDensity(true),dirtyFlagModulus(true),dirtyFlagVelocity(true),parametrisation(0),numRelaxationMechanisms(0){};
            
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
            
            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling()=0;
            
        protected:
            
            bool dirtyFlagInverseDensity; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagModulus; //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            bool dirtyFlagVelocity; //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            IndexType parametrisation; //!< ==0 if P/S-wave modulus parametrisation; ==1 Velocity-parametrisation
            
            IndexType PartitionedIn; //!< ==1 If Module is read from partitioned fileblock; ==0 if module is in single files
            IndexType PartitionedOut; //!< ==1 If Module is written to partitioned fileblock; ==0 if module is written to single files
            
            lama::DenseVector<ValueType> pWaveModulus; //!< Vector storing P-wave modulus.
            lama::DenseVector<ValueType> sWaveModulus; //!< Vector storing S-wave modulus.
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
            
            void calcModuleFromVelocity(lama::DenseVector<ValueType>& vecVelocity, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorModule );
            
            void calcVelocityFromModule(lama::DenseVector<ValueType>& vectorModule, lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity);
            
            /*! \brief Switch parameterization to velocity */
            virtual void switch2velocity()=0;
            /*! \brief Switch parameterization to modulus */
            virtual void switch2modulus()=0;
            
            /*! \brief Refresh module if they are dirty */
            virtual void refreshModule()=0;
            
            /*! \brief Refresh velocities if they are dirty */
            virtual void refreshVelocity()=0;
            
            IndexType getParametrisation();
            
            IndexType getPartitionedIn();
            IndexType getPartitionedOut();
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename, dmemo::DistributionPtr dist);
        };
    }
}

/*! \brief Getter method for parametrisation */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedIn(){
    return(PartitionedIn);
}

/*! \brief Getter method for parametrisation */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedOut(){
    return(PartitionedOut);
}

/*! \brief Getter method for parametrisation */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getParametrisation(){
    return(parametrisation);
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
    
    readModelparameter(vector,filename,dist);
    
    vector.redistribute(dist);
    
}


/*! \brief Write singe modelparameter to an external file
 *
 *  Write a single model to an external file block.
 \param vector Single modelparameter which will be written to filename
 \param filename Name of file in which modelparameter will be written
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    if (PartitionedOut==1){
        PartitionedInOut<ValueType> test;
        test.writeToDistributedFiles(vector,filename);
    } else if (PartitionedOut==0){
        vector.writeToFile(filename);
    } else  {
        COMMON_THROWEXCEPTION("Unexpected output option!")
    }
};

/*! \brief Read a modelparameter from file
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename, dmemo::DistributionPtr dist)
{
    if (PartitionedIn==1){
        PartitionedInOut::readFromDistributedFiles(vector,filename,dist);
    } else if (PartitionedIn==0) {
        PartitionedInOut::readFromOneFile(vector,filename,dist);
    } else  {
        COMMON_THROWEXCEPTION("Unexpected input option!")
    }
    
};


/*! \brief Allocate a single modelparameter
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Calculate a module from velocity
 *
 *  Calculates Module = pow(Velocity,2) *  Density
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vectorModule Modulus-Vector which is calculated
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcModuleFromVelocity(lama::DenseVector<ValueType>& vecVelocity, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorModule )
{
    
    vectorModule=vecDensity;
    vectorModule.scale(vecVelocity);
    vectorModule.scale(vecVelocity);
    
};

/*! \brief Calculate S and P wave modulus from velocities
 *  Acoustic:   pWaveModulus = rho * vP^2
 *  Elastic:    sWaveModulus = rho * vS^2
 \param vecV Velocity-Vector which will be used in the calculation (vP: Acoustic, vS: Elastic)
 \param vecDensity Density-Vector which will be used in the calculation
 \param vectorModule Modulus-Vector which is calculated
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param filenameDensity Location of external density-file which will be read in
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateModulus(lama::DenseVector<ValueType>& vecV, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vectorModule, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity)
{
    allocateModelparameter(vecV,ctx,dist);
    allocateModelparameter(vecDensity,ctx,dist);
    allocateModelparameter(vectorModule,ctx,dist);
    
    readModelparameter(vecV,filename,dist);
    readModelparameter(vecDensity,filenameDensity,dist);
    
    vecV.redistribute(dist);
    vecDensity.redistribute(dist);
    
    calcModuleFromVelocity(vecV,vecDensity,vectorModule);
    
};


/*! \brief Calculate velocities from a module
 *
 *  Calculates Velocity = sqrt( Modulu / Density )
 *
 \param vectorModule Modulus-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcVelocityFromModule(lama::DenseVector<ValueType>& vectorModule, lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecVelocity)
{
    /* Modulus = pow(velocity,2) * Density */
    /* Velocity = sqrt( Modulus / Density )  */
    vecVelocity=vecDensity;
    vecVelocity.invert(); /* = 1 / Density */
    vecVelocity.scale(vectorModule); /* = Modulus / Density */
    vecVelocity.sqrt(); /* = sqrt( Modulus / Density ) */
    
};


/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity(){
    if(dirtyFlagInverseDensity){
        dirtyFlagInverseDensity=false;
        inverseDensity.assign(density);
        inverseDensity.invert();
    }
    return(inverseDensity);
}

/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity(){
    dirtyFlagInverseDensity=true; // If density will be changed, the inverse has to be refreshed if it is accessed
    return(density);
}

/*! \brief Get reference to first Lame model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus(){
    
    // If the model is parameterized in modules, the velocity vector is now dirty
    if(parametrisation==0){
        dirtyFlagVelocity=true;
    }
    
    // If the model is parameterized in velocities AND the modulus is dirty, than recalculate
    if(dirtyFlagModulus && parametrisation==1){
        refreshModule();
    }
    
    return(pWaveModulus);
}

/*! \brief Get reference to second Lame Parameter sWaveModulus
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus(){
    
    // If the model is parameterized in modules, the velocity vector is now dirty
    if(parametrisation==0){
        dirtyFlagVelocity=true;
    }
    
    // If the model is parameterized in velocities AND the modulus is dirty, than recalculate
    if(dirtyFlagModulus && parametrisation==1){
        refreshModule();
    }
    
    return(sWaveModulus);
}

/*! \brief Get reference to P-wave velocity
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP(){
    
    // If the model is parameterized in velocities, the modulus vector is now dirty
    if(parametrisation==1){
        dirtyFlagModulus=true;
    }
    
    // If the model is parameterized in module AND the velocity is dirty, than recalculate
    if(dirtyFlagVelocity && parametrisation==0){
        refreshVelocity();
    }
    
    return(velocityP);
}

/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS(){
    
    // If the model is parameterized in velocities, the modulus vector is now dirty
    if(parametrisation==1){
        dirtyFlagModulus=true;
    }
    
    // If the model is parameterized in module AND the velocity is dirty, than recalculate
    if(dirtyFlagVelocity && parametrisation==0){
        refreshVelocity();
    }
    
    return(velocityS);
}

/*! \brief Get reference to tauP
 
 
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

