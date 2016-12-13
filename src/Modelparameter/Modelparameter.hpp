

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
            Modelparameter():dirtyFlagInverseDensity(true),dirtyFlagModulus(true),dirtyFlagAveraging(true),dirtyFlagVelocity(true),parametrisation(0),numRelaxationMechanisms(0){};
            
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
            virtual void write(std::string filename) const =0;
            
            virtual lama::DenseVector<ValueType> const& getDensity();
            virtual lama::DenseVector<ValueType> const& getDensity() const;
            virtual lama::DenseVector<ValueType> const& getInverseDensity();
            virtual lama::DenseVector<ValueType> const& getInverseDensity() const;
            virtual lama::DenseVector<ValueType> const& getPWaveModulus();
            virtual lama::DenseVector<ValueType> const& getPWaveModulus() const;
            virtual lama::DenseVector<ValueType> const& getSWaveModulus();
            virtual lama::DenseVector<ValueType> const& getSWaveModulus() const;
            virtual lama::DenseVector<ValueType> const& getVelocityP();
            virtual lama::DenseVector<ValueType> const& getVelocityP() const;
            virtual lama::DenseVector<ValueType> const& getVelocityS();
            virtual lama::DenseVector<ValueType> const& getVelocityS() const;
            
            virtual lama::DenseVector<ValueType> const& getTauP();
            virtual lama::DenseVector<ValueType> const& getTauP() const;
            virtual lama::DenseVector<ValueType> const& getTauS();
            virtual lama::DenseVector<ValueType> const& getTauS() const;
            
            virtual IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;
            
            /*! \brief Prepare the model parameters for modelling */
//            virtual void prepareForModelling()=0;
            virtual void prepareForModelling(Configuration::Configuration<ValueType> const& config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, dmemo::CommunicatorPtr comm)=0;
            
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageX();
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageX() const;
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageY();
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageY() const;
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageZ();
            virtual lama::DenseVector<ValueType> const& getInverseDensityAverageZ() const;
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageXY();
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageXY() const;
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageXZ();
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageXZ() const;
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageYZ();
            virtual lama::DenseVector<ValueType> const& getSWaveModulusAverageYZ() const;
            virtual lama::DenseVector<ValueType> const& getTauSAverageXY();
            virtual lama::DenseVector<ValueType> const& getTauSAverageXY() const;
            virtual lama::DenseVector<ValueType> const& getTauSAverageXZ();
            virtual lama::DenseVector<ValueType> const& getTauSAverageXZ() const;
            virtual lama::DenseVector<ValueType> const& getTauSAverageYZ();
            virtual lama::DenseVector<ValueType> const& getTauSAverageYZ() const;
            
            
        protected:
            
            bool dirtyFlagInverseDensity; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagModulus; //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            bool dirtyFlagAveraging; //!< ==true if averaged P/S-wave modulus has to be recalculated; ==false if averaged modulus is up to date
            bool dirtyFlagVelocity; //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            IndexType parametrisation; //!< ==0 if P/S-wave modulus parametrisation; ==1 Velocity-parametrisation
            
            lama::DenseVector<ValueType> pWaveModulus; //!< Vector storing P-wave modulus.
            lama::DenseVector<ValueType> sWaveModulus; //!< Vector storing S-wave modulus.
            lama::DenseVector<ValueType> density; //!< Vector storing Density.
            lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.
            
            lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.
            
            lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.
            
            lama::DenseVector<ValueType> inverseDensityAverageX; //!< Vector storing inverse averaged density in x-direction.
            lama::DenseVector<ValueType> inverseDensityAverageY; //!< Vector storing inverse averaged density in y-direction.
            lama::DenseVector<ValueType> inverseDensityAverageZ; //!< Vector storing inverse averaged density in z-direction.
            
            lama::DenseVector<ValueType> sWaveModulusAverageXY;  //!< Vector storing averaged s-wave modulus in xy-plan.
            lama::DenseVector<ValueType> sWaveModulusAverageXZ;  //!< Vector storing averaged s-wave modulus in xz-plan.
            lama::DenseVector<ValueType> sWaveModulusAverageYZ;  //!< Vector storing averaged s-wave modulus in yz-plan.
            
            lama::DenseVector<ValueType> tauSAverageXY;  //!< Vector storing averaged s-wave modulus in xy-plan.
            lama::DenseVector<ValueType> tauSAverageXZ;  //!< Vector storing averaged s-wave modulus in xz-plan.
            lama::DenseVector<ValueType> tauSAverageYZ;  //!< Vector storing averaged s-wave modulus in yz-plan.
            
            
            IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency; //!< Relaxation Frequency
            
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value);
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void writeModelparameter(lama::DenseVector<ValueType>const& vector, std::string filename) const;
            
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
            
            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging()=0;
            
            IndexType getParametrisation();
            
            //! \brief Initializsation of the aneraging matrices
            /*!
             *
             \param dist Distribution of the wavefield
             \param ctx Context
             \param NX Total number of grid points in X
             \param NY Total number of grid points in Y
             \param NZ Total number of grid points in Z
             \param DH Grid spacing (equidistant)
             \param DT Temporal sampling interval
             \param spatialFDorder FD-order of spatial stencils
             \param comm Communicator
             */
            virtual void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )=0;
//            
//            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
            
            void calcDensityAverageMatrixX(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixY(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            void calcSWaveModulusAverageMatrixXY(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixXZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixYZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            lama::CSRSparseMatrix<ValueType> DensityAverageMatrixX; //!< Averaging density matrix in x-direction
            lama::CSRSparseMatrix<ValueType> DensityAverageMatrixY; //!< Averaging density matrix in x-direction
            lama::CSRSparseMatrix<ValueType> DensityAverageMatrixZ; //!< Averaging density matrix in x-direction
            
            lama::CSRSparseMatrix<ValueType> sWaveModulusAverageMatrixXY; //!< Average S-wave Modulus in xy-plane
            lama::CSRSparseMatrix<ValueType> sWaveModulusAverageMatrixXZ; //!< Average S-wave Modulus in xz-plane
            lama::CSRSparseMatrix<ValueType> sWaveModulusAverageMatrixYZ; //!< Average S-wave Modulus in yz-plane
            
            void calculateInverseAveragedDensity(lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecInverseAvDensity, lama::CSRSparseMatrix<ValueType>& avDensityMatrix);
            void calculateAveragedSWaveModulus(lama::DenseVector<ValueType>& vecSWaveModulus, lama::DenseVector<ValueType>& vecAvSWaveModulus, lama::CSRSparseMatrix<ValueType>& avSWaveModulusMatrix);
            void calculateAveragedTauS(lama::DenseVector<ValueType>& vecTauS, lama::DenseVector<ValueType>& vecAvTauS, lama::CSRSparseMatrix<ValueType>& avTauSMatrix);
            
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
            
            typedef void (Modelparameter<ValueType>::*setRowElements_AvPtr)(IndexType , IndexType& , IndexType& , hmemo::WriteAccess<IndexType>& , hmemo::WriteAccess<IndexType>& ,hmemo::WriteAccess<ValueType>& , IndexType , IndexType , IndexType ); //!< Pointer to set elements functions
            
            typedef IndexType (Modelparameter<ValueType>::*calcNumberRowElements_AvPtr)(IndexType , IndexType , IndexType , IndexType); //!< Pointer to counting elements functions
            
            //! \brief Getter method for averaging density matrix in x-direction
            virtual lama::CSRSparseMatrix<ValueType>& getDensityAverageMatrixX();
            //! \brief Getter method for averaging density matrix in y-direction
            virtual lama::CSRSparseMatrix<ValueType>& getDensityAverageMatrixY();
            //! \brief Getter method for averaging density matrix in z-direction
            virtual lama::CSRSparseMatrix<ValueType>& getDensityAverageMatrixZ();
            
            //! \brief Getter method for averaging S-Wave modulus xy-plane
            virtual lama::CSRSparseMatrix<ValueType>& getSWaveModulusAverageMatrixXY();
            //! \brief Getter method for averaging S-Wave modulus xz-plane
            virtual lama::CSRSparseMatrix<ValueType>& getSWaveModulusAverageMatrixXZ();
            //! \brief Getter method for averaging S-Wave modulus yz-plane
            virtual lama::CSRSparseMatrix<ValueType>& getSWaveModulusAverageMatrixYZ();

            
            void calcAverageMatrix(lama::CSRSparseMatrix<ValueType>& Av, calcNumberRowElements_AvPtr calcNumberRowElements,setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            IndexType calcNumberRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            
            void setRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            
            void setRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            
        };
    }
}

/*! \brief Getter method for parametrisation */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getParametrisation(){
    return(parametrisation);
}

/*! \brief Getter method for relaxation frequency */
template<typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getRelaxationFrequency() const
{
    return(relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getNumRelaxationMechanisms() const
{
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(lama::DenseVector<ValueType>const& vector, std::string filename) const
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
    
    readModelparameter(vecV,filename);
    readModelparameter(vecDensity,filenameDensity);
    
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
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity(){
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
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity() const
{
//    SCAI_ASSERT(dirtyFlagInverseDensity == false, "Inverse density has to be recalculated! ");
    return(inverseDensity);
}

/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity(){
    dirtyFlagInverseDensity=true; // If density will be changed, the inverse has to be refreshed if it is accessed
    return(density);
}

/*! \brief Get reference to density model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity() const
{
    SCAI_ASSERT(dirtyFlagInverseDensity == true, "Density has to be recalculated! ");
    return(density);
}

/*! \brief Get reference to first Lame model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus(){
    
    // If the model is parameterized in modules, the velocity vector is now dirty
    if(parametrisation==0){
        dirtyFlagVelocity=true;
    }
    
    // If the model is parameterized in velocities AND the modulus is dirty, than recalculate
    if(dirtyFlagModulus && parametrisation==1){
        dirtyFlagModulus=false;
        refreshModule();
    }
    
    return(pWaveModulus);
}

/*! \brief Get reference to first Lame model parameter
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus() const
{
    SCAI_ASSERT( (dirtyFlagModulus == false) || (parametrisation == 0) , "Module has to be recalculated! ");
    return(pWaveModulus);
}

/*! \brief Get reference to second Lame Parameter sWaveModulus
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus(){
    
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

/*! \brief Get reference to second Lame Parameter sWaveModulus
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT((dirtyFlagModulus == false) || (parametrisation == 0), "Module has to be recalculated! ");
    return(sWaveModulus);
}

/*! \brief Get reference to P-wave velocity
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP(){

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

/*! \brief Get reference to P-wave velocity
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP() const
{
    SCAI_ASSERT((dirtyFlagVelocity == false) || (parametrisation == 1), "Velocity has to be recalculated! ");
    return(velocityP);
}

/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS(){
    
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

/*! \brief Get reference to S-wave velocity
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS() const
{
    SCAI_ASSERT((dirtyFlagVelocity == false) || (parametrisation == 1), "Velocity has to be recalculated! ");
    return(velocityS);
}


/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauP(){
    return(tauP);
}

/*! \brief Get reference to tauP
 *
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauP() const
{
    return(tauP);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauS(){
    return(tauS);
}

/*! \brief Get reference to tauS
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauS() const
{
    return(tauS);
}


//! \brief Function to set elements of a single row of x-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/){
    
    IndexType RowNumber_plusOne = rowNumber + 1;
    
    for(IndexType j=0; j<=1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for last element in every NX x NX matrix. Here: insert 1
        if(j==0){
            csrJALocal[countJA]=rowNumber;
            if (RowNumber_plusOne % NX == 0) {
                csrvaluesLocal[countJA]=1.0;
            } else {
                csrvaluesLocal[countJA]=1.0/2.0;
            }
            countJA++;
        } else {
            // Set elaments right to diagonal to 1.0/2.0 exept matrix elements with condition (row%NX)!=0)
            if (RowNumber_plusOne % NX != 0) {
                csrJALocal[countJA]=rowNumber+1;
                csrvaluesLocal[countJA]=1.0/2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}

//! \brief Function to set elements of a single row of y-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/){
    
    IndexType NXNY = NX * NY;
    
    
    for(IndexType j=0; j<=1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for matrix elements in the last NX x NX matrix in every submatrix NXNY x NXNY. Here: insert 1
        if(j==0){
            csrJALocal[countJA]=rowNumber;
            if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                csrvaluesLocal[countJA]=1.0;
            } else {
                csrvaluesLocal[countJA]=1.0/2.0;
            }
            countJA++;
        } else {
            // Set elaments NX right to diagonal to 1.0/2.0 not if (row%NXNY) element of {0,...,Nx-1}
            if ((rowNumber % NXNY) <= (NX * (NY - 1)-1)) {
                csrJALocal[countJA]=rowNumber+NX;
                csrvaluesLocal[countJA]=1.0/2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}


//! \brief Function to set elements of a single row of x-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
    IndexType NXNY = NX * NY;
    
    
    for(IndexType j=0; j<=1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for matrix elements in the last NXNY x NXNY matrix. Here: insert 1
        if(j==0){
            csrJALocal[countJA]=rowNumber;
            if (((rowNumber) >= (NXNY * (NZ - 1)))) {
                csrvaluesLocal[countJA]=1.0;
            } else {
                csrvaluesLocal[countJA]=1.0/2.0;
            }
            countJA++;
        } else {
            // Set elaments NXNY right to diagonal to 1.0/2.0
            if (rowNumber <= NXNY*(NZ-1)-1) {
                csrJALocal[countJA]=rowNumber+NXNY;
                csrvaluesLocal[countJA]=1.0/2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}


//! \brief Function to set elements of a single row averaging s-wave modulus matrix in xy-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/){
    
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=3; j++){
        if (j <= 1) {
            if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                // Set diagonal elements in last submatrix
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0;
                    }
                } else {
                    csrJALocal[countJA]=rowNumber + j;
                    csrvaluesLocal[countJA]=1.0/2.0;
                }
            } else {
                // Set diagonal elements and elements right next do diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    csrJALocal[countJA]=rowNumber + j;
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        } else {
            if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                // Set elements NX and NX+1 right to diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA]=rowNumber + NX;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA]=rowNumber + NX + j - 2;
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        }
        countJA++;
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}


//! \brief Function to set elements of a single row averaging s-wave modulus matrix in xz-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=3; j++){
        if (j <= 1) {
            if (rowNumber >= (NXNY * (NZ - 1))) {
                // Set diagonal elements and elements next to diagonal in last submatrix
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0;
                    }
                } else {
                    csrJALocal[countJA]=rowNumber + j;
                    csrvaluesLocal[countJA]=1.0/2.0;
                }
            } else {
                // Set diagonal elements and elements next to diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    csrJALocal[countJA]=rowNumber + j;
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                // Set elements NXNY right to diagonal
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j==2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA]=rowNumber + NXNY;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA]=rowNumber + NXNY + j-2;
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        }
        countJA++;
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}


//! \brief Function to set elements of a single row averaging s-wave modulus matrix in YZ-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
    IndexType NXNY = NX * NY;
    IndexType NXNYNZ = NX * NY * NZ;
    
    for (IndexType j=0; j<=3; j++){
        if (j <= 1) {
            if (rowNumber >= (NXNY * (NZ - 1))) {
                // Set elements in last submatrix NX x NX
                if (rowNumber >= NXNYNZ - NX) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0;
                    }
                } else {
                    // Set elements in last submatrix NXNY x NXNY
                    csrJALocal[countJA]=rowNumber + j * NX;
                    csrvaluesLocal[countJA]=1.0/2.0;
                }
            } else {
                // Set elements in diagonal submatrices NXNY x NXNY
                if ((rowNumber % NXNY) >= (NX*(NY - 1))) {
                    if (j==0) {
                        csrJALocal[countJA]=rowNumber;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    csrJALocal[countJA]=rowNumber + (j * NX);
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                // Set elements NXNY and NYNX+NX right to diagonal elements
                if((rowNumber % NXNY) >= (NX * (NY - 1))){
                    if (j==2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA]=rowNumber + NXNY;
                        csrvaluesLocal[countJA]=1.0/2.0;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA]=rowNumber + NXNY + (j-2) * NX;
                    csrvaluesLocal[countJA]=1.0/4.0;
                }
            }
        }
        countJA++;
    }
    csrIALocal[countIA]=countJA;
    countIA++;
    
}


//! \brief Calculate of averaging derivative matrix
/*!
 *
 \param Av Averaging Matrix Av
 \param calcNumberRowElements Member-Function to calculate number of elements in a single row
 \param setRowElements Member-Function to set the elements in a single row
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrix(lama::CSRSparseMatrix<ValueType>& Av, calcNumberRowElements_AvPtr calcNumberRowElements,setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist){
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); //here the local indices of each process are retrieved and stored in localIndices
    
    /* Number of grid points */
    IndexType N=NX*NY*NZ;
    
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    IndexType numLocalValues=0; // Number of local values of Matrix Df
    
    /* Calculate the number of values in each matrix */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp; // Temporary storage of the local index for the ongoing iterations
    
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        
        /* Check for elements of Av */
        numLocalValues+=(this->*calcNumberRowElements)(read_localIndices_temp,NX,NY,NZ);
        
    }
    
    /* Allocate local part to create local CSR storage*/
    hmemo::HArray<ValueType> valuesLocal(numLocalValues);
    hmemo::HArray<IndexType> csrJALocal(numLocalValues);
    hmemo::HArray<IndexType> csrIALocal(numLocalIndices+1);
    
    /* Get WriteAccess to local part */
    hmemo::WriteAccess<IndexType> write_csrJALocal(csrJALocal);
    hmemo::WriteAccess<IndexType> write_csrIALocal(csrIALocal);
    hmemo::WriteAccess<ValueType> write_valuesLocal(valuesLocal);
    
    /* Set some counters to create the CSR Storage */
    IndexType countJA=0;
    IndexType countIA=0;
    write_csrIALocal[0]=0;
    countIA++;
    
    /* Set the values into the indice arrays and the value array for CSR matrix build */
    for(IndexType i=0; i<numLocalIndices; i++){
        
        read_localIndices_temp=read_localIndices[i];
        
        // write AveragingMatrix
        (this->*setRowElements)(read_localIndices_temp,countJA,countIA,write_csrJALocal,write_csrIALocal,write_valuesLocal,NX,NY,NZ);
        
    }
    
    /* Release all read and write access */
    read_localIndices.release();
    
    write_csrJALocal.release();
    write_csrIALocal.release();
    write_valuesLocal.release();
    
    /* Create local CSR storage of Matrix D, than create distributed CSR matrix D */
    lama::CSRStorage<ValueType> Av_LocalCSR(numLocalIndices,N,numLocalValues,csrIALocal,csrJALocal,valuesLocal);
    Av_LocalCSR.compress();
    Av.assign(Av_LocalCSR,dist,dist);
    
}


//! \brief Calculate number of row elements for density averaging matrix in x-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType NX,IndexType /*NY*/, IndexType /*NZ*/){
    
    IndexType counter=0;
    
    for (IndexType j=0; j<=1; j++){
        if((rowNumber % NX) == (NX-1)){
            if (j==0){
                counter++;
            }
        } else {
            counter++;
        }
    }
    return(counter);
}


//! \brief Calculate number of row elements for density averaging matrix in y-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType NX,IndexType NY, IndexType /*NZ*/){
    
    IndexType counter=0;
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=1; j++){
        if((rowNumber % NXNY) >= (NX * (NY - 1))){
            if (j==0){
                counter++;
            }
        } else {
            counter++;
        }
    }
    return(counter);
}


//! \brief Calculate number of row elements for density averaging matrix in z-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
    IndexType counter=0;
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=1; j++){
        if(rowNumber >= (NXNY * (NZ - 1))){
            if (j==0){
                counter++;
            }
        } else {
            counter++;
        }
    }
    return(counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in xy-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType NX,IndexType NY, IndexType /*NZ*/){
    
    IndexType counter=0;
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=3; j++){
        if (j<=1) {
            if (j==0) {
                counter++;
            } else {
                if((rowNumber % NX) <= ((NX - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                if (j==2){
                    counter++;
                } else {
                    if ((rowNumber % NX) != (NX - 1)) {
                        counter++;
                    }
                }
            }
        }
    }
    return(counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in xz-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
    IndexType counter=0;
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=3; j++){
        if (j<=1) {
            if (j==0) {
                counter++;
            } else {
                if((rowNumber % NX) <= ((NX - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                if (j==2){
                    counter++;
                } else {
                    if ((rowNumber % NX) != (NX - 1)) {
                        counter++;
                    }
                }
            }
            
        }
    }
    return(counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in yz-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template<typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
    IndexType counter=0;
    IndexType NXNY = NX * NY;
    
    for (IndexType j=0; j<=3; j++){
        if (j<=1) {
            if (j==0) {
                counter++;
            } else {
                if((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                if (j==2){
                    counter++;
                } else {
                    if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                        counter++;
                    }
                }
            }
            
        }
    }
    return(counter);
}



//! \brief Calculate density averaging matrix in x-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixX(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixX, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixX, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixX, NX, NY, NZ, dist);
}


//! \brief Calculate density averaging matrix in y-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixY(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixY, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixY, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixY, NX, NY, NZ, dist);
}


//! \brief Calculate density averaging matrix in z-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixZ, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixZ, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixZ, NX, NY, NZ, dist);
}


//! \brief Calculate s-wave modulus averaging matrix in x-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXY(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixXY, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXY, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXY, NX, NY, NZ, dist);
}


//! \brief Calculate s-wave modulus averaging matrix in y-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ, NX, NY, NZ, dist);
}


//! \brief Calculate s-wave modulus averaging matrix in z-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixYZ(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixYZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixYZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixYZ, NX, NY, NZ, dist);
}


/*! \brief calculate averaged s-wave modulus
 *
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateInverseAveragedDensity(lama::DenseVector<ValueType>& vecDensity, lama::DenseVector<ValueType>& vecInverseAvDensity, lama::CSRSparseMatrix<ValueType>& avDensityMatrix)
{
    vecInverseAvDensity = avDensityMatrix* vecDensity;
    vecInverseAvDensity.invert();
}


/*! \brief calculate averaged s-wave modulus
 *
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedSWaveModulus(lama::DenseVector<ValueType>& vecSWaveModulus, lama::DenseVector<ValueType>& vecAvSWaveModulus, lama::CSRSparseMatrix<ValueType>& avSWaveModulusMatrix)
{
    vecAvSWaveModulus = vecSWaveModulus;
    vecAvSWaveModulus.invert();
    vecAvSWaveModulus = avSWaveModulusMatrix * vecAvSWaveModulus;
    vecAvSWaveModulus.invert();
}


/*! \brief calculate averaged tauS
 *
 \param vecAvTauS Averaged tauS vector which is calculated
 \param avTauSMatrix Averaging matrix which is used to calculate averaged vector
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedTauS(lama::DenseVector<ValueType>& vecTauS, lama::DenseVector<ValueType>& vecAvTauS, lama::CSRSparseMatrix<ValueType>& avTauSMatrix)
{
    vecAvTauS = vecTauS;
    vecAvTauS = avTauSMatrix * vecAvTauS;
}


//! \brief Getter method for averaging density matrix in x-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityAverageMatrixX(){
    return(DensityAverageMatrixX);
}


//! \brief Getter method for averaging density matrix in y-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityAverageMatrixY(){
    return(DensityAverageMatrixY);
}


//! \brief Getter method for averaging density matrix in z-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityAverageMatrixZ(){
    return(DensityAverageMatrixZ);
}


//! \brief Getter method for averaging S-wave modulus matrix x-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageMatrixXY(){
    return(sWaveModulusAverageMatrixXY);
}

//! \brief Getter method for averaging S-wave modulus matrix y-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageMatrixXZ(){
    return(sWaveModulusAverageMatrixXZ);
}

//! \brief Getter method for averaging S-wave modulus matrix z-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageMatrixYZ(){
    return(sWaveModulusAverageMatrixYZ);
}

/*! \brief Get reference to averaged density in x-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(inverseDensityAverageX);
}

/*! \brief Get reference to averaged density in y-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(inverseDensityAverageY);
}

/*! \brief Get reference to averaged density in z-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(inverseDensityAverageZ);
}

/*! \brief Get reference to averaged s-wave modulus in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(sWaveModulusAverageXY);
}

/*! \brief Get reference to averaged s-wave modulus in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(sWaveModulusAverageXZ);
}

/*! \brief Get reference to averaged s-wave modulus in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(sWaveModulusAverageYZ);
}

/*! \brief Get reference to averaged tauS in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(tauSAverageXY);
}

/*! \brief Get reference to averaged tauS in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(tauSAverageXZ);
}

/*! \brief Get reference to averaged tauS in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ(){
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if(dirtyFlagAveraging == true){
        calculateAveraging();
    }
    return(tauSAverageYZ);
}

/*! \brief Get reference to averaged density in x-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(inverseDensityAverageX);
}

/*! \brief Get reference to averaged density in y-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(inverseDensityAverageY);
}

/*! \brief Get reference to averaged density in z-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(inverseDensityAverageZ);
}

/*! \brief Get reference to averaged s-wave modulus in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(sWaveModulusAverageXY);
}

/*! \brief Get reference to averaged s-wave modulus in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(sWaveModulusAverageXZ);
}

/*! \brief Get reference to averaged s-wave modulus in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(sWaveModulusAverageYZ);
}

/*! \brief Get reference to averaged tauS in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(tauSAverageXY);
}

/*! \brief Get reference to averaged tauS in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(tauSAverageXZ);
}

/*! \brief Get reference to averaged tauS in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>const& KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return(tauSAverageYZ);
}


