

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
            
            //! \brief Getter method for averaging density matrix in x-direction
            virtual lama::CSRSparseMatrix<ValueType>& get_xAvDensityMatrix();
            //! \brief Getter method for averaging density matrix in y-direction
            virtual lama::CSRSparseMatrix<ValueType>& get_yAvDensityMatrix();
            //! \brief Getter method for averaging density matrix in z-direction
            virtual lama::CSRSparseMatrix<ValueType>& get_zAvDensityMatrix();
            
            //! \brief Getter method for averaging S-Wave modulus xy-plane
            virtual lama::CSRSparseMatrix<ValueType>& get_xyAvSWaveModulusMatrix();
            //! \brief Getter method for averaging S-Wave modulus xz-plane
            virtual lama::CSRSparseMatrix<ValueType>& get_xzAvSWaveModulusMatrix();
            //! \brief Getter method for averaging S-Wave modulus yz-plane
            virtual lama::CSRSparseMatrix<ValueType>& get_yzAvSWaveModulusMatrix();
            
            virtual lama::DenseVector<ValueType>& get_InverseXAvDensity();
            virtual lama::DenseVector<ValueType>& get_InverseYAvDensity();
            virtual lama::DenseVector<ValueType>& get_InverseZAvDensity();
            virtual lama::DenseVector<ValueType>& get_xyAvSWaveModulus();
            virtual lama::DenseVector<ValueType>& get_xzAvSWaveModulus();
            virtual lama::DenseVector<ValueType>& get_yzAvSWaveModulus();
            virtual lama::DenseVector<ValueType>& get_xyAvTauS();
            virtual lama::DenseVector<ValueType>& get_xzAvTauS();
            virtual lama::DenseVector<ValueType>& get_yzAvTauS();
            
            
            
            
        protected:
            
            bool dirtyFlagInverseDensity; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagModulus; //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
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
            
            lama::DenseVector<ValueType> inverseXAvDensity; //!< Vector storing inverse averaged density in x-direction.
            lama::DenseVector<ValueType> inverseYAvDensity; //!< Vector storing inverse averaged density in y-direction.
            lama::DenseVector<ValueType> inverseZAvDensity; //!< Vector storing inverse averaged density in z-direction.
            
            lama::DenseVector<ValueType> xyAvSWaveModulus;  //!< Vector storing averaged s-wave modulus in xy-plan.
            lama::DenseVector<ValueType> xzAvSWaveModulus;  //!< Vector storing averaged s-wave modulus in xz-plan.
            lama::DenseVector<ValueType> yzAvSWaveModulus;  //!< Vector storing averaged s-wave modulus in yz-plan.
            
            lama::DenseVector<ValueType> xyAvTauS;  //!< Vector storing averaged s-wave modulus in xy-plan.
            lama::DenseVector<ValueType> xzAvTauS;  //!< Vector storing averaged s-wave modulus in xz-plan.
            lama::DenseVector<ValueType> yzAvTauS;  //!< Vector storing averaged s-wave modulus in yz-plan.
            
            
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
            
            void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
            
            void calc_xAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calc_yAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calc_zAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            void calc_xyAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calc_xzAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            void calc_yzAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            lama::CSRSparseMatrix<ValueType> xAvDensityMatrix; //!< Averaging density matrix in x-direction
            lama::CSRSparseMatrix<ValueType> yAvDensityMatrix; //!< Averaging density matrix in x-direction
            lama::CSRSparseMatrix<ValueType> zAvDensityMatrix; //!< Averaging density matrix in x-direction
            
            lama::CSRSparseMatrix<ValueType> xyAvSWaveModulusMatrix; //!< Average S-wave Modulus in xy-plane
            lama::CSRSparseMatrix<ValueType> xzAvSWaveModulusMatrix; //!< Average S-wave Modulus in xz-plane
            lama::CSRSparseMatrix<ValueType> yzAvSWaveModulusMatrix; //!< Average S-wave Modulus in yz-plane
            
            void calculateInverseAveragedDensity(lama::DenseVector<ValueType>& vecInverseAvDensity, lama::CSRSparseMatrix<ValueType>& avDensityMatrix);
            void calculateAveragedSWaveModulus(lama::DenseVector<ValueType>& vecAvSWaveModulus, lama::CSRSparseMatrix<ValueType>& avSWaveModulusMatrix);
            void calculateAveragedTauS(lama::DenseVector<ValueType>& vecAvTauS, lama::CSRSparseMatrix<ValueType>& avTauSMatrix);
            
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
            
            typedef void (Modelparameter<ValueType>::*setRowElements_AvPtr)(IndexType , IndexType& , IndexType& , hmemo::WriteAccess<IndexType>& , hmemo::WriteAccess<IndexType>& ,hmemo::WriteAccess<ValueType>& , IndexType , IndexType , IndexType ); //!< Pointer to set elements functions
            
            typedef IndexType (Modelparameter<ValueType>::*calcNumberRowElements_AvPtr)(IndexType , IndexType , IndexType , IndexType); //!< Pointer to counting elements functions
            
            void calc_AvMatrix(lama::CSRSparseMatrix<ValueType>& Av, calcNumberRowElements_AvPtr calcNumberRowElements,setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist);
            
            IndexType calcNumberRowElements_xAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_yAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_zAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_xyAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_xzAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_yzAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ);
            
            void setRowElements_xAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_yAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_zAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            
            void setRowElements_xyAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_xzAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_yzAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_xAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_yAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_zAvDensityMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_xyAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_xzAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_yzAvSWaveModulusMatrix(IndexType rowNumber, IndexType& countJA, IndexType& countIA, hmemo::WriteAccess<IndexType>& csrJALocal, hmemo::WriteAccess<IndexType>& csrIALocal,hmemo::WriteAccess<ValueType>& csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_AvMatrix(lama::CSRSparseMatrix<ValueType>& Av, calcNumberRowElements_AvPtr calcNumberRowElements,setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist){
    
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



//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    initializeMatrices(dist,ctx,config.getNX(), config.getNY(), config.getNZ(), config.getDH(), config.getDT(), comm);
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_xAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType /*NY*/, IndexType /*NZ*/){
    
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_yAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType /*NZ*/){
    
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_zAvDensityMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_xyAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType /*NZ*/){
    
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_xzAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
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
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_yzAvSWaveModulusMatrix(IndexType rowNumber, IndexType NX,IndexType NY, IndexType NZ){
    
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_xAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(xAvDensityMatrix, &Modelparameter<ValueType>::calcNumberRowElements_xAvDensityMatrix, &Modelparameter<ValueType>::setRowElements_xAvDensityMatrix, NX, NY, NZ, dist);
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_yAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(yAvDensityMatrix, &Modelparameter<ValueType>::calcNumberRowElements_yAvDensityMatrix, &Modelparameter<ValueType>::setRowElements_yAvDensityMatrix, NX, NY, NZ, dist);
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_zAvDensityMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(zAvDensityMatrix, &Modelparameter<ValueType>::calcNumberRowElements_zAvDensityMatrix, &Modelparameter<ValueType>::setRowElements_zAvDensityMatrix, NX, NY, NZ, dist);
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_xyAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(xyAvSWaveModulusMatrix, &Modelparameter<ValueType>::calcNumberRowElements_xyAvSWaveModulusMatrix, &Modelparameter<ValueType>::setRowElements_xyAvSWaveModulusMatrix, NX, NY, NZ, dist);
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_xzAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(xzAvSWaveModulusMatrix, &Modelparameter<ValueType>::calcNumberRowElements_xzAvSWaveModulusMatrix, &Modelparameter<ValueType>::setRowElements_xzAvSWaveModulusMatrix, NX, NY, NZ, dist);
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc_yzAvSWaveModulusMatrix(IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist)
{
    calc_AvMatrix(yzAvSWaveModulusMatrix, &Modelparameter<ValueType>::calcNumberRowElements_yzAvSWaveModulusMatrix, &Modelparameter<ValueType>::setRowElements_yzAvSWaveModulusMatrix, NX, NY, NZ, dist);
}


/*! \brief calculate averaged s-wave modulus
 *
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateInverseAveragedDensity(lama::DenseVector<ValueType>& vecInverseAvDensity, lama::CSRSparseMatrix<ValueType>& avDensityMatrix)
{
    vecInverseAvDensity = avDensityMatrix* density;
    vecInverseAvDensity.invert();
}


/*! \brief calculate averaged s-wave modulus
 *
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedSWaveModulus(lama::DenseVector<ValueType>& vecAvSWaveModulus, lama::CSRSparseMatrix<ValueType>& avSWaveModulusMatrix)
{
    vecAvSWaveModulus = sWaveModulus;
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedTauS(lama::DenseVector<ValueType>& vecAvTauS, lama::CSRSparseMatrix<ValueType>& avTauSMatrix)
{
    vecAvTauS = tauS;
    vecAvTauS = avTauSMatrix * vecAvTauS;
}


//! \brief Getter method for averaging density matrix in x-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xAvDensityMatrix(){
    return(xAvDensityMatrix);
}


//! \brief Getter method for averaging density matrix in y-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_yAvDensityMatrix(){
    return(yAvDensityMatrix);
}


//! \brief Getter method for averaging density matrix in z-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_zAvDensityMatrix(){
    return(zAvDensityMatrix);
}


//! \brief Getter method for averaging S-wave modulus matrix x-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xyAvSWaveModulusMatrix(){
    return(xyAvSWaveModulusMatrix);
}

//! \brief Getter method for averaging S-wave modulus matrix y-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xzAvSWaveModulusMatrix(){
    return(xzAvSWaveModulusMatrix);
}

//! \brief Getter method for averaging S-wave modulus matrix z-direction
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_yzAvSWaveModulusMatrix(){
    return(yzAvSWaveModulusMatrix);
}

/*! \brief Get reference to averaged density in x-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_InverseXAvDensity(){
    return(inverseXAvDensity);
}

/*! \brief Get reference to averaged density in y-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_InverseYAvDensity(){
    return(inverseYAvDensity);
}

/*! \brief Get reference to averaged density in z-direction
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_InverseZAvDensity(){
    return(inverseZAvDensity);
}

/*! \brief Get reference to averaged s-wave modulus in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xyAvSWaveModulus(){
    return(xyAvSWaveModulus);
}

/*! \brief Get reference to averaged s-wave modulus in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xzAvSWaveModulus(){
    return(xzAvSWaveModulus);
}

/*! \brief Get reference to averaged s-wave modulus in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_yzAvSWaveModulus(){
    return(yzAvSWaveModulus);
}

/*! \brief Get reference to averaged tauS in xy-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xyAvTauS(){
    return(xyAvTauS);
}

/*! \brief Get reference to averaged tauS in xz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_xzAvTauS(){
    return(xzAvTauS);
}

/*! \brief Get reference to averaged tauS in yz-plane
 */
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Modelparameter::Modelparameter<ValueType>::get_yzAvTauS(){
    return(yzAvTauS);
}

