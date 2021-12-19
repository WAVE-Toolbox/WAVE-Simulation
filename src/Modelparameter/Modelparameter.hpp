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
            Modelparameter(){};

            //! Default destructor.
            ~Modelparameter(){};

            //! \brief Modelparameter pointer
            typedef std::shared_ptr<Modelparameter<ValueType>> ModelparameterPtr;

            /* Common */
            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryModel(scai::dmemo::DistributionPtr dist);

            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) = 0;

            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void write(std::string filename, scai::IndexType fileFormat) const = 0;
            virtual void writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) = 0;

            ValueType const getCenterFrequencyCPML() const;
            virtual std::string getEquationType() const = 0;
            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual std::vector<ValueType> getRelaxationFrequency() const;
            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;
            virtual scai::lama::Vector<ValueType> const &getReflectivity() const;
            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation);            
            virtual void setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity);
            
            virtual void resetReflectivity();
            virtual void calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT) = 0;

            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms) = 0;
            virtual void setRelaxationFrequency(std::vector<ValueType> const setRelaxationFrequency) = 0;
            
            virtual void calcRockMatrixParameter(Configuration::Configuration const &config) = 0;             
            virtual void calcPetrophysicsFromWaveModulus() = 0;   
            virtual void calcWaveModulusFromPetrophysics() = 0; 
            
            scai::IndexType getParameterisation() const;
            void setParameterisation(scai::IndexType const setParameterisation);
            bool getEffectiveParameterisation() const;
            void setEffectiveParameterisation(bool const setEffectiveParameterisation);
            scai::IndexType getInversionType() const;
            void setInversionType(scai::IndexType const setInversionType);
            scai::IndexType getGradientType() const;
            void setGradientType(scai::IndexType const setGradientType);
            scai::IndexType getDecomposeType() const;
            void setDecomposeType(scai::IndexType const setDecomposeType);
            
            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;
            virtual void prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) = 0;

            virtual void applyThresholds(Configuration::Configuration const &config) = 0;
            
            virtual void getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) = 0;

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate);            
            scai::lama::SparseVector<ValueType> getShrinkVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidthLeft, scai::IndexType boundaryWidthRight);
            
            virtual void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;

            /* Operator overloading */
            KITGPI::Modelparameter::Modelparameter<ValueType> &operator=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);
            KITGPI::Modelparameter::Modelparameter<ValueType> &operator-=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);
            KITGPI::Modelparameter::Modelparameter<ValueType> &operator+=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs);

            /* Seismic */
            virtual scai::lama::Vector<ValueType> const &getDensity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() = 0;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() const = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() const = 0;
            virtual scai::lama::Vector<ValueType> const &getDensityRockMatrix() const = 0;
            virtual scai::lama::Vector<ValueType> const &getBulkModulusRockMatrix() const = 0;
            virtual scai::lama::Vector<ValueType> const &getShearModulusRockMatrix() const = 0;

            virtual scai::lama::Vector<ValueType> const &getTauP() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauS() const = 0;
            
            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity) = 0;
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP) = 0;
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS) = 0;
            virtual void setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix) = 0;
            virtual void setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix) = 0;
            virtual void setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix) = 0;
            
            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP) = 0;
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS) = 0;
            
            virtual scai::lama::DenseVector<ValueType> const &getBiotCoefficient() = 0;
            virtual scai::lama::DenseVector<ValueType> const &getBiotCoefficient() const = 0; 
            virtual scai::lama::DenseVector<ValueType> const &getBulkModulusM() = 0;
            virtual scai::lama::DenseVector<ValueType> const &getBulkModulusM() const = 0;
            virtual scai::lama::DenseVector<ValueType> const &getBulkModulusKf() = 0;
            virtual scai::lama::DenseVector<ValueType> const &getBulkModulusKf() const = 0;
            
            virtual ValueType const getDensityWater() const = 0;
            virtual ValueType const getDensityAir() const = 0;
            virtual ValueType const getBulkModulusWater() const = 0; 
            virtual ValueType const getBulkModulusAir() const = 0;       
            virtual ValueType const getCriticalPorosity() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX() const = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ() const = 0;
                        
            /* EM */
            virtual scai::lama::Vector<ValueType> const &getMagneticPermeability() const = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() = 0;
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() const = 0;
            
            virtual scai::lama::DenseVector<ValueType> const getElectricConductivityRealEffective() const = 0;
            virtual scai::lama::DenseVector<ValueType> const getDielectricPermittivityRealEffective() const = 0;
            virtual scai::lama::DenseVector<ValueType> const getElectricConductivityStatic(scai::lama::DenseVector<ValueType> const dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> const electricConductivityRealEffective, scai::IndexType calculateType) = 0;
            virtual scai::lama::DenseVector<ValueType> const getDielectricPermittivityStatic(scai::lama::DenseVector<ValueType> const dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> const electricConductivityRealEffective, scai::IndexType calculateType) = 0;
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const = 0;

            virtual scai::lama::Vector<ValueType> const &getElectricConductivityWater() const = 0;
            virtual scai::lama::Vector<ValueType> const &getRelativeDieletricPeimittivityRockMatrix() const = 0;    
            
            virtual void setMagneticPermeability(scai::lama::Vector<ValueType> const &setMagneticPermeability) = 0;
            virtual void setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity) = 0;
            virtual void setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity) = 0;
            
            virtual void setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity) = 0;
            virtual void setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity) = 0;
            
            virtual void setElectricConductivityWater(scai::lama::Vector<ValueType> const &setElectricConductivityWater) = 0;
            virtual void setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix) = 0; 
            
            virtual ValueType const getDielectricPermittivityVacuum() const = 0; 
            virtual ValueType const getElectricConductivityReference() const = 0; 
            virtual ValueType const getTauDielectricPermittivityReference() const = 0; 
            virtual ValueType const getTauElectricConductivityReference() const = 0; 
            virtual ValueType const getRelativeDielectricPermittivityWater() const = 0; 
            virtual ValueType const getRelativeDielectricPermittivityVacuum() const = 0; 
            
            virtual void setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const) = 0;  
            virtual ValueType const getArchie_a() const = 0; 
            virtual ValueType const getArchie_m() const = 0;  
            virtual ValueType const getArchie_n() const = 0;            
            virtual void calcElectricConductivityReference(ValueType const CenterFrequencyCPML) = 0;
            
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ() const = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY() = 0;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX() = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX() const = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY() = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX() = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX() const = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY() = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageX() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageX() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageY() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageZ() const = 0;
            
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY() const = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ() = 0;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ() const = 0;

          protected:
            /* Common */
            void calcInverseAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
            void calcAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
            
            virtual void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecModulus) = 0;
            virtual void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecModulus, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecVelocity) = 0;

            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

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
            
            bool dirtyFlagAveraging = true;      //!< ==true if averaged P/S-wave modulus has to be recalculated; ==false if averaged modulus is up to date
            
            scai::IndexType parameterisation = 0;
            bool effectiveParameterisation = 0;
            scai::IndexType inversionType = 0;
            scai::IndexType gradientType = 0;
            scai::IndexType decomposeType = 0;
            scai::IndexType fileFormat;      //!< 1=mtx 2=lmf
            std::string equationType;
            
            scai::lama::DenseVector<ValueType> porosity; 
            scai::lama::DenseVector<ValueType> saturation; 
            scai::lama::DenseVector<ValueType> reflectivity; //!< Vector storing reflectivity.
            
            scai::IndexType numRelaxationMechanisms = 0; //!< Number of relaxation mechanisms
            std::vector<ValueType> relaxationFrequency;           //!< Relaxation Frequency
            ValueType centerFrequencyCPML;           //!< Relaxation Frequency

            SparseFormat averageMatrixX;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixY;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixZ;                                 //!< Averaging density matrix in x-direction
            SparseFormat averageMatrixXY;                                //!< Average S-wave Modulus in xy-plane
            SparseFormat averageMatrixXZ;                                //!< Average S-wave Modulus in xz-plane
            SparseFormat averageMatrixYZ;                                //!< Average S-wave Modulus in yz-plane

            /* Seismic */
            virtual void calcAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix) = 0;
            
            bool dirtyFlagInverseDensity = true; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            bool dirtyFlagPWaveModulus = true;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            bool dirtyFlagSWaveModulus = true;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            
            scai::lama::DenseVector<ValueType> pWaveModulus;   //!< Vector storing P-wave modulus.
            scai::lama::DenseVector<ValueType> sWaveModulus;   //!< Vector storing S-wave modulus.
            scai::lama::DenseVector<ValueType> density;        //!< Vector storing Density.
            scai::lama::DenseVector<ValueType> inverseDensity; //!< Vector storing inverted density.

            scai::lama::DenseVector<ValueType> velocityP; //!< Vector storing P-wave velocity.
            scai::lama::DenseVector<ValueType> velocityS; //!< Vector storing S-wave velocity.

            scai::lama::DenseVector<ValueType> tauP; //!< Vector storing tauP for viscoelastic modelling.
            scai::lama::DenseVector<ValueType> tauS; //!< Vector storing tauS for viscoelastic modelling.

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

            /* EM */
            bool dirtyFlagVelocivityEM = true;   //!< ==true if EM-wave velocity has to be recalculated; ==false if velocity is up to date
            
            ValueType aArchie;
            ValueType mArchie;
            ValueType nArchie;   
            ValueType const MagneticPermeabilityVacuum = 1.2566370614e-6;       
            ValueType const DielectricPermittivityVacuum = 8.8541878176e-12;   
            ValueType ElectricConductivityReference;   
            ValueType const TauDielectricPermittivityReference = 0.01; 
            ValueType const RelativeDielectricPermittivityWater = 81; 
            ValueType const RelativeDielectricPermittivityVacuum = 1;
            
            scai::lama::DenseVector<ValueType> velocivityEM;   //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> magneticPermeability;        //!< Vector storing MagneticPermeability.
            scai::lama::DenseVector<ValueType> electricConductivity; //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> dielectricPermittivity; //!< Vector storing dielectricPermittivity.
            
            scai::lama::DenseVector<ValueType> tauElectricConductivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            
            scai::lama::DenseVector<ValueType> electricConductivityWater; //!< Vector electricConductivityWater.
            scai::lama::DenseVector<ValueType> relativeDieletricPeimittivityRockMatrix; //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityAverageYZ; //!< Vector storing inverse averaged magneticPermeability in yz-plan.
            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityAverageXZ; //!< Vector storing inverse averaged magneticPermeability in xz-plan.
            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityAverageXY; //!< Vector storing inverse averaged magneticPermeability in xy-plan.
            scai::lama::DenseVector<ValueType> electricConductivityAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> electricConductivityAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> electricConductivityAverageZ; //!< Vector storing averaged modulus in z-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityAverageZ; //!< Vector storing averaged modulus in z-direction.
            
            scai::lama::DenseVector<ValueType> tauElectricConductivityAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> tauElectricConductivityAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> tauElectricConductivityAverageZ; //!< Vector storing averaged modulus in z-direction.              
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageZ; //!< Vector storing averaged modulus in z-direction.

            std::vector<ValueType> relaxationTime;     // relaxation time of electric displacement
            
          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
        };
    }
}
