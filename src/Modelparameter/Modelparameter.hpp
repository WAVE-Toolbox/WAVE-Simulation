

#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "../Common/HostPrint.hpp"
#include "../Configuration/Configuration.hpp"
#include "../PartitionedInOut/PartitionedInOut.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! \brief Abstract class for a single Modelparameter (Subsurface properties)
        /*!
         * This class handels a single modelparameter.
         * As this class is an abstract class, all constructors are protected.
         */
        template <typename ValueType>
        class Modelparameter
        {
          public:
            //! Default constructor.
            Modelparameter() : dirtyFlagInverseDensity(true), dirtyFlagAveraging(true), dirtyFlagPWaveModulus(true),dirtyFlagSWaveModulus(true),parametrisation(0), numRelaxationMechanisms(0){};

            //! Default destructor.
            ~Modelparameter(){};

            //! \brief Modelparameter pointer
            typedef std::shared_ptr<Modelparameter<ValueType>> ModelparameterPtr;

            /*! \brief Abstract initialization function
             * Standard initialisation function
             \param ctx Context
             \param dist Distribution
             \param filename filename to read modelparameters (endings will be added by derived classes)
             \param partitionedIn Partitioned input
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) = 0;

            /*! \brief Abstract initialisation function
             * Standard initialisation function
             \param config Configuration from configuration file
             \param ctx Context
             \param dist Distribution
             */
            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write modelparameters (endings will be added by derived classes)
             \param partitionedOut Partitioned output
             */
            virtual void write(std::string filename, IndexType partitionedOut) const = 0;

            virtual scai::lama::Vector const &getDensity() const;
            virtual scai::lama::Vector const &getInverseDensity();
            virtual scai::lama::Vector const &getInverseDensity() const;
            virtual scai::lama::Vector const &getPWaveModulus();
            virtual scai::lama::Vector const &getPWaveModulus() const;
            virtual scai::lama::Vector const &getSWaveModulus();
            virtual scai::lama::Vector const &getSWaveModulus() const;
            virtual scai::lama::Vector const &getVelocityP() const;
            virtual scai::lama::Vector const &getVelocityS() const;

            virtual scai::lama::Vector const &getTauP() const;
            virtual scai::lama::Vector const &getTauS() const;

            virtual IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setDensity(scai::lama::Vector const &setDensity);
            virtual void setVelocityP(scai::lama::Vector const &setVelocityP);
            virtual void setVelocityS(scai::lama::Vector const &setVelocityS);

            virtual void setTauP(scai::lama::Vector const &setTauP);
            virtual void setTauS(scai::lama::Vector const &setTauS);

            virtual void setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;

            virtual scai::lama::Vector const &getInverseDensityAverageX();
            virtual scai::lama::Vector const &getInverseDensityAverageX() const;
            virtual scai::lama::Vector const &getInverseDensityAverageY();
            virtual scai::lama::Vector const &getInverseDensityAverageY() const;
            virtual scai::lama::Vector const &getInverseDensityAverageZ();
            virtual scai::lama::Vector const &getInverseDensityAverageZ() const;
            virtual scai::lama::Vector const &getSWaveModulusAverageXY();
            virtual scai::lama::Vector const &getSWaveModulusAverageXY() const;
            virtual scai::lama::Vector const &getSWaveModulusAverageXZ();
            virtual scai::lama::Vector const &getSWaveModulusAverageXZ() const;
            virtual scai::lama::Vector const &getSWaveModulusAverageYZ();
            virtual scai::lama::Vector const &getSWaveModulusAverageYZ() const;
            virtual scai::lama::Vector const &getTauSAverageXY();
            virtual scai::lama::Vector const &getTauSAverageXY() const;
            virtual scai::lama::Vector const &getTauSAverageXZ();
            virtual scai::lama::Vector const &getTauSAverageXZ() const;
            virtual scai::lama::Vector const &getTauSAverageYZ();
            virtual scai::lama::Vector const &getTauSAverageYZ() const;
	    
	    virtual bool getDirtyFlagPWaveModulus() const;
	    virtual bool getDirtyFlagSWaveModulus() const;
	    virtual bool getDirtyFlagInverseDensity() const;
	    virtual bool getDirtyFlagAveraging() const;

	    virtual void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;

	    /* Operator overloading */
            KITGPI::Modelparameter::Modelparameter<ValueType> &operator=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);
	    KITGPI::Modelparameter::Modelparameter<ValueType> &operator-=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);
	    KITGPI::Modelparameter::Modelparameter<ValueType> &operator+=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);
	    
          protected:
            bool dirtyFlagInverseDensity; 	//!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagAveraging;      	//!< ==true if averaged P/S-wave modulus has to be recalculated; ==false if averaged modulus is up to date      
            bool dirtyFlagPWaveModulus;        //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            bool dirtyFlagSWaveModulus;        //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
        
            IndexType parametrisation;    //!< ==0 if P/S-wave modulus parametrisation; ==1 Velocity-parametrisation
            IndexType PartitionedIn;  //!< ==1 If Modulus is read from partitioned fileblock; ==0 if modulus is in single files
            IndexType PartitionedOut; //!< ==1 If Modulus is written to partitioned fileblock; ==0 if modulus is written to single files

            scai::lama::DenseVector<ValueType> pWaveModulus;   //!< Vector storing P-wave modulus.
            scai::lama::DenseVector<ValueType> sWaveModulus;   //!< Vector storing S-wave modulus.
            scai::lama::DenseVector<ValueType> density;        //!< Vector storing Density.
            scai::lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for visco-elastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for visco-elastic modelling.

            scai::lama::DenseVector<ValueType> inverseDensityAverageX; //!< Vector storing inverse averaged density in x-direction.
            scai::lama::DenseVector<ValueType> inverseDensityAverageY; //!< Vector storing inverse averaged density in y-direction.
            scai::lama::DenseVector<ValueType> inverseDensityAverageZ; //!< Vector storing inverse averaged density in z-direction.

            scai::lama::DenseVector<ValueType> sWaveModulusAverageXY; //!< Vector storing averaged s-wave modulus in xy-plan.
            scai::lama::DenseVector<ValueType> sWaveModulusAverageXZ; //!< Vector storing averaged s-wave modulus in xz-plan.
            scai::lama::DenseVector<ValueType> sWaveModulusAverageYZ; //!< Vector storing averaged s-wave modulus in yz-plan.

            scai::lama::DenseVector<ValueType> tauSAverageXY; //!< Vector storing averaged s-wave modulus in xy-plan.
            scai::lama::DenseVector<ValueType> tauSAverageXZ; //!< Vector storing averaged s-wave modulus in xz-plan.
            scai::lama::DenseVector<ValueType> tauSAverageYZ; //!< Vector storing averaged s-wave modulus in yz-plan.

            IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;     //!< Relaxation Frequency

            void initModelparameter(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar value);
            void initModelparameter(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            void writeModelparameter(scai::lama::Vector const &vector, std::string filename, IndexType partitionedOut) const;

            void calcModulusFromVelocity(scai::lama::Vector &vecVelocity, scai::lama::Vector &vecDensity, scai::lama::Vector &vectorModulus);

            void calcVelocityFromModulus(scai::lama::Vector &vectorModulus, scai::lama::Vector &vecV, scai::lama::Vector &vecDensity);

            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging() = 0;

            IndexType getParametrisation();
            IndexType getPartitionedIn();
            IndexType getPartitionedOut();
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
             \param comm Communicator
             */
            virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm) = 0;

            void calcDensityAverageMatrixX(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixY(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

            void calcSWaveModulusAverageMatrixXY(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixXZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixYZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat DensityAverageMatrixX;                          //!< Averaging density matrix in x-direction
            SparseFormat DensityAverageMatrixY;                          //!< Averaging density matrix in x-direction
            SparseFormat DensityAverageMatrixZ;                          //!< Averaging density matrix in x-direction
            SparseFormat sWaveModulusAverageMatrixXY;                    //!< Average S-wave Modulus in xy-plane
            SparseFormat sWaveModulusAverageMatrixXZ;                    //!< Average S-wave Modulus in xz-plane
            SparseFormat sWaveModulusAverageMatrixYZ;                    //!< Average S-wave Modulus in yz-plane

            void calculateInverseAveragedDensity(scai::lama::Vector &vecDensity, scai::lama::Vector &vecInverseAvDensity, scai::lama::Matrix &avDensityMatrix);
            void calculateAveragedSWaveModulus(scai::lama::Vector &vecSWaveModulus, scai::lama::Vector &vecAvSWaveModulus, scai::lama::Matrix &avSWaveModulusMatrix);
            void calculateAveragedTauS(scai::lama::Vector &vecTauS, scai::lama::Vector &vecAvTauS, scai::lama::Matrix &avTauSMatrix);

          private:
            void allocateModelparameter(scai::lama::Vector &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void readModelparameter(scai::lama::Vector &vector, std::string filename, scai::dmemo::DistributionPtr dist, IndexType partitionedIn);

            typedef void (Modelparameter<ValueType>::*setRowElements_AvPtr)(IndexType, IndexType &, IndexType &, scai::hmemo::WriteAccess<IndexType> &, scai::hmemo::WriteAccess<IndexType> &, scai::hmemo::WriteAccess<ValueType> &, IndexType, IndexType, IndexType); //!< Pointer to set elements functions

            typedef IndexType (Modelparameter<ValueType>::*calcNumberRowElements_AvPtr)(IndexType, IndexType, IndexType, IndexType); //!< Pointer to counting elements functions

            void calcAverageMatrix(scai::lama::Matrix &Av, calcNumberRowElements_AvPtr calcNumberRowElements, setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

            IndexType calcNumberRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
            IndexType calcNumberRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);

            void setRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);

            void setRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
        };
    }
}
