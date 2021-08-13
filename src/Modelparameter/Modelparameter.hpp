

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
#include "../ForwardSolver/Derivatives/Derivatives.hpp"

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
            Modelparameter() : dirtyFlagInverseDensity(true), dirtyFlagAveraging(true), dirtyFlagPWaveModulus(true), dirtyFlagSWaveModulus(true), parameterisation(0), inversionType(0), numRelaxationMechanisms(0){};

            //! Default destructor.
            ~Modelparameter(){};

            //! \brief Modelparameter pointer
            typedef std::shared_ptr<Modelparameter<ValueType>> ModelparameterPtr;

            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryModel(scai::dmemo::DistributionPtr dist);

            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) = 0;

            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void write(std::string filename, scai::IndexType fileFormat) const = 0;
            virtual void writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) = 0;

            virtual std::string getEquationType() const = 0;

            virtual scai::lama::Vector<ValueType> const &getDensity() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity();
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() const;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() const;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus();
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() const;
            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;
            virtual scai::lama::Vector<ValueType> const &getDensityRockMatrix() const;
            virtual scai::lama::Vector<ValueType> const &getBulkModulusRockMatrix() const;
            virtual scai::lama::Vector<ValueType> const &getShearModulusRockMatrix() const;
            virtual scai::lama::Vector<ValueType> const &getReflectivity() const;

            virtual scai::lama::Vector<ValueType> const &getTauP() const;
            virtual scai::lama::Vector<ValueType> const &getTauS() const;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity);
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP);
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS);
            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation);
            virtual void setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix);
            virtual void setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix);
            virtual void setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix);
            virtual void setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity);
            virtual void resetReflectivity();
            virtual void calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT) = 0;

            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP);
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS);

            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);
            virtual void calcRockMatrixParameter(Configuration::Configuration const &config) = 0; 
            virtual void calcWaveModulusFromPetrophysics() = 0;      
            virtual void calcPetrophysicsFromWaveModulus() = 0;        
            scai::lama::DenseVector<ValueType> const &getBiotCoefficient(); 
            scai::lama::DenseVector<ValueType> const &getBiotCoefficient() const; 
            scai::lama::DenseVector<ValueType> const &getBulkModulusM();
            scai::lama::DenseVector<ValueType> const &getBulkModulusM() const;
            scai::lama::DenseVector<ValueType> const &getBulkModulusKf();
            scai::lama::DenseVector<ValueType> const &getBulkModulusKf() const;
            
            ValueType const getDensityWater() const; 
            ValueType const getDensityAir() const; 
            ValueType const getBulkModulusWater() const; 
            ValueType const getBulkModulusAir() const; 
            ValueType const getCriticalPorosity() const; 

            scai::IndexType getParameterisation() const;
            void setParameterisation(scai::IndexType const setParameterisation);
            scai::IndexType getInversionType() const;
            void setInversionType(scai::IndexType const setInversionType);
            scai::IndexType getGradientType() const;
            void setGradientType(scai::IndexType const setGradientType);
            scai::IndexType getDecomposeType() const;
            void setDecomposeType(scai::IndexType const setDecomposeType);
            
            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;
            void prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

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
            
            virtual void getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) = 0;
            virtual void setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) = 0;
            
            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate);            
            scai::lama::SparseVector<ValueType> getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth);
            
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

            scai::IndexType parameterisation;
            scai::IndexType inversionType;
            scai::IndexType gradientType;
            scai::IndexType decomposeType;
            scai::IndexType fileFormat;      //!< 1=mtx 2=lmf

            std::string equationType;

            scai::lama::DenseVector<ValueType> pWaveModulus;   //!< Vector storing P-wave modulus.
            scai::lama::DenseVector<ValueType> sWaveModulus;   //!< Vector storing S-wave modulus.
            scai::lama::DenseVector<ValueType> density;        //!< Vector storing Density.
            scai::lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for viscoelastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for viscoelastic modelling.

            scai::lama::DenseVector<ValueType> porosity; 
            scai::lama::DenseVector<ValueType> saturation; 
            scai::lama::DenseVector<ValueType> reflectivity; //!< Vector storing reflectivity.
            scai::lama::DenseVector<ValueType> densityRockMatrix;        //!< Vector storing Density.
            scai::lama::DenseVector<ValueType> bulkModulusRockMatrix;   //!< Vector storing P-wave modulus.
            scai::lama::DenseVector<ValueType> shearModulusRockMatrix;   //!< Vector storing S-wave modulus.

            scai::lama::DenseVector<ValueType> bulkModulusKf; 
            scai::lama::DenseVector<ValueType> BiotCoefficient;
            scai::lama::DenseVector<ValueType> bulkModulusM;
            
            ValueType const DensityWater = 1000;           
            ValueType const DensityAir = 1.29; 
            ValueType const VelocityPWater = 1473;           
            ValueType const VelocityPAir = 340; 
            ValueType const CriticalPorosity = 0.4; 
            
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

            void calcVelocityFromModulus(scai::lama::DenseVector<ValueType> &vectorModulus, scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecVelocity);

            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging() = 0;

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

            void calcInverseAveragedParameter(scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecInverseAvDensity, scai::lama::Matrix<ValueType> &avDensityMatrix);
            void calcAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix);
            void calcAveragedParameter(scai::lama::Vector<ValueType> &vecTauS, scai::lama::Vector<ValueType> &vecAvTauS, scai::lama::Matrix<ValueType> &avTauSMatrix);

          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
        };
    }
}
