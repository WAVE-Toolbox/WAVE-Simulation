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
        class ModelparameterEM
        {
          public:
            //! Default constructor.
            ModelparameterEM() : dirtyFlagAveraging(true), dirtyFlagVelocivityEM(true), dirtyFlagElectricConductivityOptical(true), dirtyFlagDielectricPermittivityOptical(true), parameterisation(0), inversionType(0), numRelaxationMechanisms(0){};

            //! Default destructor.
            ~ModelparameterEM(){};

            //! \brief Modelparameter pointer
            typedef std::shared_ptr<ModelparameterEM<ValueType>> ModelparameterPtr;

            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryModel(scai::dmemo::DistributionPtr dist);

            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) = 0;

            virtual void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void write(std::string filename, scai::IndexType fileFormat) const = 0;
            void writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) const;

            virtual std::string getEquationType() const = 0;

            virtual scai::lama::Vector<ValueType> const &getMagneticPermeability() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivity() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivity() const;
            
            virtual scai::lama::Vector<ValueType> const &getVelocityEM();
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() const;
            
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOptical();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOptical() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOptical();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOptical() const;
            virtual ValueType const &getTauElectricDisplacement();
            virtual ValueType const &getTauElectricDisplacement() const;
            virtual scai::lama::Vector<ValueType> const &getTauElectricConductivity() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const;

            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityWater() const;
            virtual scai::lama::Vector<ValueType> const &getRelativeDieletricPeimittivityRockMatrix() const;
            virtual scai::lama::Vector<ValueType> const &getReflectivity() const;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setMagneticPermeability(scai::lama::Vector<ValueType> const &setMagneticPermeability);
            virtual void setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity);
            virtual void setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity);
            
            virtual void setElectricConductivityOptical(scai::lama::Vector<ValueType> const &setElectricConductivityOptical);
            virtual void setDielectricPermittivityOptical(scai::lama::Vector<ValueType> const &setDielectricPermittivityOptical);
            virtual void setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity);
            virtual void setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity);
            
            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation);
            virtual void setElectricConductivityWater(scai::lama::Vector<ValueType> const &setElectricConductivityWater);
            virtual void setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix); 
            virtual void setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity);
            virtual void resetReflectivity();
            virtual void calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT);
            
            ValueType const getDielectricPermittivityVacuum() const; 
            ValueType const getElectricConductivityReference() const; 
            ValueType const getTauDielectricPermittivityReference() const; 
            ValueType const getTauElectricConductivityReference() const; 
            ValueType const getRelativeDielectricPermittivityWater() const; 
            ValueType const getRelativeDielectricPermittivityVacuum() const; 
            
            void setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const);  
            ValueType const getArchie_a() const; 
            ValueType const getArchie_m() const;  
            ValueType const getArchie_n() const;            
            void calcRockMatrixParameter(Configuration::Configuration const &config);
            void calcElectricConductivityReference(ValueType const CenterFrequencyCPML);
            void calcWaveModulusFromPetrophysics();
            void calcPetrophysicsFromWaveModulus();
            
            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);
            
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

            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageYZ() const;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXZ() const;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityAverageXY() const;
            
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageX();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageY();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageZ();
            virtual scai::lama::Vector<ValueType> const &getElectricConductivityOpticalAverageZ() const;            
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageX();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageY();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageZ();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityOpticalAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityAverageZ() const;
            
            virtual void getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) = 0;
            virtual void setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) = 0;

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Declare Sparse-Matrix
            SparseFormat getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate);            
            scai::lama::SparseVector<ValueType> getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth);
            
            virtual void minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) = 0;
            virtual void plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) = 0;
            virtual void assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) = 0;

            /* Operator overloading */
            KITGPI::Modelparameter::ModelparameterEM<ValueType> &operator=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ModelparameterEM<ValueType> &operator-=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ModelparameterEM<ValueType> &operator+=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs);

          protected:
            bool dirtyFlagAveraging;      //!< ==true if averaged EM-wave modulus has to be recalculated; ==false if averaged modulus is up to date
            bool dirtyFlagVelocivityEM;   //!< ==true if EM-wave velocity has to be recalculated; ==false if velocity is up to date
            bool dirtyFlagElectricConductivityOptical; //!< ==true if electricConductivityOptical has to be recalulated;
            bool dirtyFlagDielectricPermittivityOptical; //!< ==true if dielectricPermittivityOptical has to be recalulated;

            scai::IndexType parameterisation;
            scai::IndexType inversionType;
            scai::IndexType gradientType;
            scai::IndexType decomposeType;
            scai::IndexType fileFormat;      //!< 1=mtx 2=lmf

            ValueType aArchie;
            ValueType mArchie;
            ValueType nArchie;   
            ValueType const MagneticPermeabilityVacuum = 1.2566370614e-6;       
            ValueType const DielectricPermittivityVacuum = 8.8541878176e-12;   
            ValueType ElectricConductivityReference;   
            ValueType const TauDielectricPermittivityReference = 0.01; 
            ValueType const RelativeDielectricPermittivityWater = 81; 
            ValueType const RelativeDielectricPermittivityVacuum = 1;
            
            std::string equationType;

            scai::lama::DenseVector<ValueType> velocivityEM;   //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> magneticPermeability;        //!< Vector storing MagneticPermeability.
            scai::lama::DenseVector<ValueType> electricConductivity; //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> dielectricPermittivity; //!< Vector storing dielectricPermittivity.
            
            scai::lama::DenseVector<ValueType> electricConductivityOptical; //!< Vector storing dielectricPermittivityOptical for visco-emem modelling.
            scai::lama::DenseVector<ValueType> dielectricPermittivityOptical; //!< Vector storing dielectricPermittivityOptical for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauElectricConductivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivity; //!< Vector storing tauDielectricPermittivity for visco-emem modelling.
            
            scai::lama::DenseVector<ValueType> porosity; //!< Vector storing porosity.
            scai::lama::DenseVector<ValueType> saturation; //!< Vector storing saturation.
            scai::lama::DenseVector<ValueType> reflectivity; //!< Vector storing reflectivity.
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
            
            scai::lama::DenseVector<ValueType> electricConductivityOpticalAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> electricConductivityOpticalAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> electricConductivityOpticalAverageZ; //!< Vector storing averaged modulus in z-direction.          
            scai::lama::DenseVector<ValueType> dielectricPermittivityOpticalAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityOpticalAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityOpticalAverageZ; //!< Vector storing averaged modulus in z-direction.              
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityAverageZ; //!< Vector storing averaged modulus in z-direction.
         
            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;           //!< Relaxation Frequency
            ValueType tauElectricDisplacement;     // relaxation time of electric displacement

            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecDielectricPermittivity);

            void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecDielectricPermittivity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecVelocity);

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

            SparseFormat averageMatrixX;                                 //!< Averaging magneticPermeability matrix in x-direction
            SparseFormat averageMatrixY;                                 //!< Averaging magneticPermeability matrix in y-direction
            SparseFormat averageMatrixZ;                                 //!< Averaging magneticPermeability matrix in z-direction
            SparseFormat averageMatrixXY;                                //!< Average EM-wave Modulus in xy-plane
            SparseFormat averageMatrixXZ;                                //!< Average EM-wave Modulus in xz-plane
            SparseFormat averageMatrixYZ;                                //!< Average EM-wave Modulus in yz-plane

            void calcInverseAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
            void calcAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
                         
          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
        };
    }
}
