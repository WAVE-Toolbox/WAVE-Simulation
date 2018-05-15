

#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
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
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType partitionedIn) = 0;

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
            virtual void write(std::string filename, scai::IndexType partitionedOut) const = 0;

            virtual scai::lama::Vector<ValueType> const &getDensity() const;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity();
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() const;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const;

            virtual scai::lama::Vector<ValueType> const &getTauP() const;
            virtual scai::lama::Vector<ValueType> const &getTauS() const;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity);
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP);
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS);

            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP);
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS);

            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);

            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;

            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX();
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY();
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ();
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() const;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY();
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY() const;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ();
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ() const;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ();
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ() const;
	    
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
        
            scai::IndexType parametrisation;    //!< ==0 if P/S-wave modulus parametrisation; ==1 Velocity-parametrisation
            scai::IndexType PartitionedIn;  //!< ==1 If Modulus is read from partitioned fileblock; ==0 if modulus is in single files
            scai::IndexType PartitionedOut; //!< ==1 If Modulus is written to partitioned fileblock; ==0 if modulus is written to single files

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

            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;     //!< Relaxation Frequency

            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType partitionedIn);

            void writeModelparameter(scai::lama::Vector<ValueType> const &vector, std::string filename, scai::IndexType partitionedOut) const;

            void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vectorModulus);

            void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vectorModulus, scai::lama::Vector<ValueType> &vecV, scai::lama::Vector<ValueType> &vecDensity);

            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging() = 0;

            scai::IndexType getParametrisation();
            scai::IndexType getPartitionedIn();
            scai::IndexType getPartitionedOut();
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
            virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm) = 0;

            void calcDensityAverageMatrixX(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixY(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcDensityAverageMatrixZ(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

            void calcSWaveModulusAverageMatrixXY(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixXZ(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
            void calcSWaveModulusAverageMatrixYZ(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat DensityAverageMatrixX;                          //!< Averaging density matrix in x-direction
            SparseFormat DensityAverageMatrixY;                          //!< Averaging density matrix in x-direction
            SparseFormat DensityAverageMatrixZ;                          //!< Averaging density matrix in x-direction
            SparseFormat sWaveModulusAverageMatrixXY;                    //!< Average S-wave Modulus in xy-plane
            SparseFormat sWaveModulusAverageMatrixXZ;                    //!< Average S-wave Modulus in xz-plane
            SparseFormat sWaveModulusAverageMatrixYZ;                    //!< Average S-wave Modulus in yz-plane

            void calculateInverseAveragedDensity(scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecInverseAvDensity, scai::lama::Matrix<ValueType> &avDensityMatrix);
            void calculateAveragedSWaveModulus(scai::lama::Vector<ValueType> &vecSWaveModulus, scai::lama::Vector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix);
            void calculateAveragedTauS(scai::lama::Vector<ValueType> &vecTauS, scai::lama::Vector<ValueType> &vecAvTauS, scai::lama::Matrix<ValueType> &avTauSMatrix);

          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void readModelparameter(scai::lama::Vector<ValueType> &vector, std::string filename, scai::dmemo::DistributionPtr dist, scai::IndexType partitionedIn);

            typedef void (Modelparameter<ValueType>::*setRowElements_AvPtr)(scai::IndexType, scai::IndexType &, scai::IndexType &, scai::hmemo::WriteAccess<scai::IndexType> &, scai::hmemo::WriteAccess<scai::IndexType> &, scai::hmemo::WriteAccess<ValueType> &, scai::IndexType, scai::IndexType, scai::IndexType); //!< Pointer to set elements functions

            typedef scai::IndexType (Modelparameter<ValueType>::*calcNumberRowElements_AvPtr)(scai::IndexType, scai::IndexType, scai::IndexType, scai::IndexType); //!< Pointer to counting elements functions

            void calcAverageMatrix(scai::lama::Matrix<ValueType> &Av, calcNumberRowElements_AvPtr calcNumberRowElements, setRowElements_AvPtr setRowElements, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

            scai::IndexType calcNumberRowElements_DensityAverageMatrixX(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType calcNumberRowElements_DensityAverageMatrixY(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType calcNumberRowElements_DensityAverageMatrixZ(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType calcNumberRowElements_SWaveModulusAverageMatrixXY(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType calcNumberRowElements_SWaveModulusAverageMatrixXZ(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType calcNumberRowElements_SWaveModulusAverageMatrixYZ(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            void setRowElements_DensityAverageMatrixX(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            void setRowElements_DensityAverageMatrixY(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            void setRowElements_DensityAverageMatrixZ(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            void setRowElements_SWaveModulusAverageMatrixXY(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixXZ(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            void setRowElements_SWaveModulusAverageMatrixYZ(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
        };
    }
}
