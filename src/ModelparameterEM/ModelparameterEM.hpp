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
#include "../Modelparameter/Modelparameter.hpp"

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
        class ModelparameterEM : public Modelparameter<ValueType>
        {
          public:
            //! Default constructor.
            ModelparameterEM(){};

            //! Default destructor.
            ~ModelparameterEM(){};

            /* Common */
            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;
            
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) = 0;

            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void write(std::string filename, scai::IndexType fileFormat) const = 0;
            void writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) override;

            virtual std::string getEquationType() const = 0;

            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;
            void prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) override;

            virtual void applyThresholds(Configuration::Configuration const &config) = 0;

            void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms) override;
            void setRelaxationFrequency(std::vector<ValueType> const setRelaxationFrequency) override;
            virtual void calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT) override;
            
            virtual void calcRockMatrixParameter(Configuration::Configuration const &config) override; 
            virtual void calcWaveModulusFromPetrophysics() override;      
            virtual void calcPetrophysicsFromWaveModulus() override;   
            
            virtual void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) = 0;

            /* Seismic */
            virtual scai::lama::Vector<ValueType> const &getDensity() const override;
            virtual scai::lama::Vector<ValueType> const &getVelocityP() const override;
            virtual scai::lama::Vector<ValueType> const &getVelocityS() const override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensity() const override;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() override;
            virtual scai::lama::Vector<ValueType> const &getPWaveModulus() const override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulus() const override;
            virtual scai::lama::Vector<ValueType> const &getDensityRockMatrix() const override;
            virtual scai::lama::Vector<ValueType> const &getBulkModulusRockMatrix() const override;
            virtual scai::lama::Vector<ValueType> const &getShearModulusRockMatrix() const override;

            virtual scai::lama::Vector<ValueType> const &getTauP() const override;
            virtual scai::lama::Vector<ValueType> const &getTauS() const override;
            
            virtual void setDensity(scai::lama::Vector<ValueType> const &setDensity) override;
            virtual void setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP) override;
            virtual void setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS) override;
            virtual void setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix) override;
            virtual void setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix) override;
            virtual void setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix) override;
            
            virtual void setTauP(scai::lama::Vector<ValueType> const &setTauP) override;
            virtual void setTauS(scai::lama::Vector<ValueType> const &setTauS) override;
     
            scai::lama::DenseVector<ValueType> const &getBiotCoefficient() override; 
            scai::lama::DenseVector<ValueType> const &getBiotCoefficient() const override; 
            scai::lama::DenseVector<ValueType> const &getBulkModulusM() override;
            scai::lama::DenseVector<ValueType> const &getBulkModulusM() const override;
            scai::lama::DenseVector<ValueType> const &getBulkModulusKf() override;
            scai::lama::DenseVector<ValueType> const &getBulkModulusKf() const override;
            
            ValueType const getDensityWater() const override; 
            ValueType const getDensityAir() const override; 
            ValueType const getBulkModulusWater() const override; 
            ValueType const getBulkModulusAir() const override; 
            ValueType const getCriticalPorosity() const override; 

            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX() override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageX() const override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY() override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageY() const override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ() override;
            virtual scai::lama::Vector<ValueType> const &getInverseDensityAverageZ() const override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() const override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() const override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() override;
            virtual scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() const override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY() override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXY() const override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ() override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageXZ() const override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ() override;
            virtual scai::lama::Vector<ValueType> const &getTauSAverageYZ() const override;
            
            /* EM */
            virtual scai::lama::Vector<ValueType> const &getMagneticPermeability() const override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() const override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() const override;
            
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() override;
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() const override;
            
            virtual scai::lama::DenseVector<ValueType> const getElectricConductivityRealEffective() const;
            virtual scai::lama::DenseVector<ValueType> const getDielectricPermittivityRealEffective() const;
            virtual scai::lama::DenseVector<ValueType> const getElectricConductivityStatic(scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> electricConductivityRealEffective);
            virtual scai::lama::DenseVector<ValueType> const getDielectricPermittivityStatic(scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> electricConductivityRealEffective);
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() const override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const override;

            virtual scai::lama::Vector<ValueType> const &getElectricConductivityWater() const override;
            virtual scai::lama::Vector<ValueType> const &getRelativeDieletricPeimittivityRockMatrix() const override;

            virtual void setMagneticPermeability(scai::lama::Vector<ValueType> const &setMagneticPermeability) override;
            virtual void setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity) override;
            virtual void setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity) override;
            
            virtual void setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity) override;
            virtual void setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity) override;
            
            virtual void setElectricConductivityWater(scai::lama::Vector<ValueType> const &setElectricConductivityWater) override;
            virtual void setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix) override; 
            
            ValueType const getDielectricPermittivityVacuum() const override; 
            ValueType const getElectricConductivityReference() const override; 
            ValueType const getTauDielectricPermittivityReference() const override; 
            ValueType const getTauElectricConductivityReference() const override; 
            ValueType const getRelativeDielectricPermittivityWater() const override; 
            ValueType const getRelativeDielectricPermittivityVacuum() const override; 
            
            void setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const) override;  
            ValueType const getArchie_a() const override; 
            ValueType const getArchie_m() const override;  
            ValueType const getArchie_n() const override;            
            void calcElectricConductivityReference(ValueType const CenterFrequencyCPML) override;
                        
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ() override;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ() const override;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ() override;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ() const override;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY() override;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY() const override;
            
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX() override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX() const override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY() override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY() const override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ() override;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ() const override;
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX() override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX() const override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY() override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY() const override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ() override;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ() const override;
            
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageX() override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageX() const override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageY() override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageY() const override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageZ() override;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivityAverageZ() const override;
            
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX() override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX() const override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY() override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY() const override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ() override;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ() const override;
            
          protected:              
            /* Common */
            void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecDielectricPermittivity) override;

            void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecDielectricPermittivity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecVelocity) override;

            /*! \brief Calculate Averaging if they are required */
            virtual void calculateAveraging() = 0;
            
            virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) = 0;
            virtual void purgeMatrices() = 0;

            using Modelparameter<ValueType>::dirtyFlagAveraging;      //!< ==true if averaged EM-wave modulus has to be recalculated; ==false if averaged modulus is up to date
            
            using Modelparameter<ValueType>::porosity; //!< Vector storing porosity.
            using Modelparameter<ValueType>::saturation; //!< Vector storing saturation.
            using Modelparameter<ValueType>::reflectivity; //!< Vector storing reflectivity.
            
            using Modelparameter<ValueType>::parameterisation;
            using Modelparameter<ValueType>::inversionType;
            using Modelparameter<ValueType>::gradientType;
            using Modelparameter<ValueType>::decomposeType;
            using Modelparameter<ValueType>::fileFormat;      //!< 1=mtx 2=lmf

            using Modelparameter<ValueType>::equationType;
            using Modelparameter<ValueType>::numRelaxationMechanisms; //!< Number of relaxation mechanisms
            using Modelparameter<ValueType>::relaxationFrequency;           //!< Relaxation Frequency
            using Modelparameter<ValueType>::centerFrequencyCPML;
            
            using Modelparameter<ValueType>::averageMatrixX;                                 //!< Averaging magneticPermeability matrix in x-direction
            using Modelparameter<ValueType>::averageMatrixY;                                 //!< Averaging magneticPermeability matrix in y-direction
            using Modelparameter<ValueType>::averageMatrixZ;                                 //!< Averaging magneticPermeability matrix in z-direction
            using Modelparameter<ValueType>::averageMatrixXY;                                //!< Average EM-wave Modulus in xy-plane
            using Modelparameter<ValueType>::averageMatrixXZ;                                //!< Average EM-wave Modulus in xz-plane
            using Modelparameter<ValueType>::averageMatrixYZ;                                //!< Average EM-wave Modulus in yz-plane
            
            /* Seismic */
            void calcAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix) override;
            using Modelparameter<ValueType>::dirtyFlagInverseDensity; //!< ==true if inverseDensity has to be recalulated; ==false if inverseDensity is up to date
            using Modelparameter<ValueType>::dirtyFlagPWaveModulus;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            using Modelparameter<ValueType>::dirtyFlagSWaveModulus;   //!< ==true if P/S-wave modulus has to be recalculated; ==false if modulus is up to date
            
            using Modelparameter<ValueType>::pWaveModulus;   //!< Vector storing P-wave modulus.
            using Modelparameter<ValueType>::sWaveModulus;   //!< Vector storing S-wave modulus.
            using Modelparameter<ValueType>::density;        //!< Vector storing Density.
            using Modelparameter<ValueType>::inverseDensity; //!< Vector storing inverted density.

            using Modelparameter<ValueType>::velocityP; //!< Vector storing P-wave velocity.
            using Modelparameter<ValueType>::velocityS; //!< Vector storing S-wave velocity.

            using Modelparameter<ValueType>::tauP; //!< Vector storing tauP for viscoelastic modelling.
            using Modelparameter<ValueType>::tauS; //!< Vector storing tauS for viscoelastic modelling.

            using Modelparameter<ValueType>::densityRockMatrix;        //!< Vector storing Density.
            using Modelparameter<ValueType>::bulkModulusRockMatrix;   //!< Vector storing P-wave modulus.
            using Modelparameter<ValueType>::shearModulusRockMatrix;   //!< Vector storing S-wave modulus.

            using Modelparameter<ValueType>::bulkModulusKf; 
            using Modelparameter<ValueType>::BiotCoefficient;
            using Modelparameter<ValueType>::bulkModulusM;
            
            using Modelparameter<ValueType>::DensityWater;           
            using Modelparameter<ValueType>::DensityAir; 
            using Modelparameter<ValueType>::VelocityPWater;           
            using Modelparameter<ValueType>::VelocityPAir; 
            using Modelparameter<ValueType>::CriticalPorosity; 
            
            using Modelparameter<ValueType>::inverseDensityAverageX; //!< Vector storing inverse averaged density in x-direction.
            using Modelparameter<ValueType>::inverseDensityAverageY; //!< Vector storing inverse averaged density in y-direction.
            using Modelparameter<ValueType>::inverseDensityAverageZ; //!< Vector storing inverse averaged density in z-direction.

            using Modelparameter<ValueType>::sWaveModulusAverageXY; //!< Vector storing averaged s-wave modulus in xy-plan.
            using Modelparameter<ValueType>::sWaveModulusAverageXZ; //!< Vector storing averaged s-wave modulus in xz-plan.
            using Modelparameter<ValueType>::sWaveModulusAverageYZ; //!< Vector storing averaged s-wave modulus in yz-plan.

            using Modelparameter<ValueType>::tauSAverageXY; //!< Vector storing averaged s-wave modulus in xy-plan.
            using Modelparameter<ValueType>::tauSAverageXZ; //!< Vector storing averaged s-wave modulus in xz-plan.
            using Modelparameter<ValueType>::tauSAverageYZ; //!< Vector storing averaged s-wave modulus in yz-plan.
                        
            /* EM */
            using Modelparameter<ValueType>::dirtyFlagVelocivityEM;   //!< ==true if EM-wave velocity has to be recalculated; ==false if velocity is up to date

            using Modelparameter<ValueType>::aArchie;
            using Modelparameter<ValueType>::mArchie;
            using Modelparameter<ValueType>::nArchie;   
            using Modelparameter<ValueType>::MagneticPermeabilityVacuum;       
            using Modelparameter<ValueType>::DielectricPermittivityVacuum;   
            using Modelparameter<ValueType>::ElectricConductivityReference;   
            using Modelparameter<ValueType>::TauDielectricPermittivityReference; 
            using Modelparameter<ValueType>::RelativeDielectricPermittivityWater; 
            using Modelparameter<ValueType>::RelativeDielectricPermittivityVacuum;
            
            using Modelparameter<ValueType>::velocivityEM;   //!< Vector storing EM-wave velocity.
            using Modelparameter<ValueType>::magneticPermeability;        //!< Vector storing MagneticPermeability.
            using Modelparameter<ValueType>::electricConductivity; //!< Vector storing EM-wave velocity.
            using Modelparameter<ValueType>::dielectricPermittivity; //!< Vector storing dielectricPermittivity.
            
            using Modelparameter<ValueType>::tauElectricConductivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            using Modelparameter<ValueType>::tauDielectricPermittivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            
            using Modelparameter<ValueType>::electricConductivityWater; //!< Vector electricConductivityWater.
            using Modelparameter<ValueType>::relativeDieletricPeimittivityRockMatrix; //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageYZ; //!< Vector storing inverse averaged magneticPermeability in yz-plan.
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXZ; //!< Vector storing inverse averaged magneticPermeability in xz-plan.
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXY; //!< Vector storing inverse averaged magneticPermeability in xy-plan.
            using Modelparameter<ValueType>::electricConductivityAverageX; //!< Vector storing averaged modulus in x-direction.
            using Modelparameter<ValueType>::electricConductivityAverageY; //!< Vector storing averaged modulus in y-direction.
            using Modelparameter<ValueType>::electricConductivityAverageZ; //!< Vector storing averaged modulus in z-direction.
            using Modelparameter<ValueType>::dielectricPermittivityAverageX; //!< Vector storing averaged modulus in x-direction.
            using Modelparameter<ValueType>::dielectricPermittivityAverageY; //!< Vector storing averaged modulus in y-direction.
            using Modelparameter<ValueType>::dielectricPermittivityAverageZ; //!< Vector storing averaged modulus in z-direction.
            
            using Modelparameter<ValueType>::tauElectricConductivityAverageX; //!< Vector storing averaged modulus in x-direction.
            using Modelparameter<ValueType>::tauElectricConductivityAverageY; //!< Vector storing averaged modulus in y-direction.
            using Modelparameter<ValueType>::tauElectricConductivityAverageZ; //!< Vector storing averaged modulus in z-direction.          
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageX; //!< Vector storing averaged modulus in x-direction.
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageY; //!< Vector storing averaged modulus in y-direction.
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageZ; //!< Vector storing averaged modulus in z-direction.
         
            using Modelparameter<ValueType>::relaxationTime;     // relaxation time of electric displacement

        };
    }
}
