

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

#include "../Acquisition/Coordinates.hpp"
#include "../Common/Common.hpp"
#include "../Common/HostPrint.hpp"
#include "../Configuration/Configuration.hpp"
#include <iostream>

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
            Modelparameter() : dirtyFlagInverseDensity(true), dirtyFlagAveraging(true), dirtyFlagPWaveModulus(true), dirtyFlagSWaveModulus(true), parametrisation(0), numRelaxationMechanisms(0){};

            //! Default destructor.
            ~Modelparameter(){};

            //! \brief Modelparameter pointer
            typedef std::shared_ptr<Modelparameter<ValueType>> ModelparameterPtr;

            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryModel(scai::dmemo::DistributionPtr dist);

            /*! \brief Abstract initialization function
             * Standard initialisation function
             \param ctx Context
             \param dist Distribution
             \param filename filename to read modelparameters (endings will be added by derived classes)
             \param fileFormat Input file format 1=mtx 2=lmf
             */
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) = 0;

            /*! \brief Abstract initialisation function
             * Standard initialisation function
             \param config Configuration from configuration file
             \param ctx Context
             \param dist Distribution
             */
            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            /*! \brief Abstract write function
             *
             * Standard write function
             *
             \param filename filename to write modelparameters (endings will be added by derived classes)
             \param fileFormat Output file format 1=mtx 2=lmf
             */
            virtual void write(std::string filename, scai::IndexType fileFormat) const = 0;

            virtual std::string getEquationType() const = 0;

            virtual scai::lama::DenseVector<ValueType> const &getDensity() const;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity();
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() const;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() const;
            virtual scai::lama::DenseVector<ValueType> const &getVelocityP() const;
            virtual scai::lama::DenseVector<ValueType> const &getVelocityS() const;

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
            virtual void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;

            virtual void applyThresholds(Configuration::Configuration const &config) = 0;

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
            
            virtual void getModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &modelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd)=0;
            virtual void setModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &invertedModelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NY, scai::IndexType NXBig, scai::IndexType NYBig, scai::IndexType boundaryWidth)=0;
            
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
            bool dirtyFlagInverseDensity; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagAveraging;      //!< ==true if averaged P/S-wave modulus has to be recalculated; ==false if averaged modulus is up to date
            bool dirtyFlagPWaveModulus;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            bool dirtyFlagSWaveModulus;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date

            scai::IndexType parametrisation; //!< ==0 if P/S-wave modulus parametrisation; ==1 Velocity-parametrisation
            scai::IndexType fileFormat;      //!< 1=mtx 2=lmf

            std::string equationType;

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
            ValueType relaxationFrequency;           //!< Relaxation Frequency

            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vectorModulus);

            void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vectorModulus, scai::lama::Vector<ValueType> &vecV, scai::lama::Vector<ValueType> &vecDensity);

            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging() = 0;

            scai::IndexType getParametrisation();

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoordinates);
            
            scai::lama::SparseVector<ValueType> getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D &cutCoordinates, scai::IndexType NX, scai::IndexType NYBig, scai::IndexType boundaryWidth);
            scai::lama::DenseVector<ValueType> smoothParameter(Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, scai::lama::DenseVector<ValueType> parameter, Acquisition::coordinate3D &cutCoordinates, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NXBig, scai::IndexType NYBig);
            
            //! \brief Initializsation of the aneraging matrices
            /*!
             *
             \param dist Distribution of the wavefield
             \param ctx Context
             \param comm Communicator
             */
            virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) = 0;
            virtual void purgeMatrices() = 0;

            void calcAverageMatrixX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
            void calcAverageMatrixY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
            void calcAverageMatrixZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);

            void calc4PointAverageMatrixRow(scai::IndexType rowIndex, scai::IndexType pX[], scai::IndexType pY[], scai::IndexType pz[], scai::lama::MatrixAssembly<ValueType> &assembly, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            void calcAverageMatrixXY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
            void calcAverageMatrixXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
            void calcAverageMatrixYZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
            
            SparseFormat averageMatrixX;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixY;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixZ;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixXY;                                //!< Average S-wave Modulus in xy-plane
            SparseFormat averageMatrixXZ;                                //!< Average S-wave Modulus in xz-plane
            SparseFormat averageMatrixYZ;                                //!< Average S-wave Modulus in yz-plane

            void calculateInverseAveragedDensity(scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecInverseAvDensity, scai::lama::Matrix<ValueType> &avDensityMatrix);
            void calculateAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix);
            void calculateAveragedTauS(scai::lama::Vector<ValueType> &vecTauS, scai::lama::Vector<ValueType> &vecAvTauS, scai::lama::Matrix<ValueType> &avTauSMatrix);

          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
        };
    }
}
