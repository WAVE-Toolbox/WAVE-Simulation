

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
        class ModelparameterEM
        {
          public:
            //! Default constructor.
            ModelparameterEM() : dirtyFlagAveraging(true), dirtyFlagVelocivityEM(true), dirtyFlagConductivityEMoptical(true), dirtyFlagDielectricPermittivityEMoptical(true), parameterisation(0), inversionType(0), numRelaxationMechanisms(0){};

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

            virtual scai::lama::Vector<ValueType> const &getMagneticPermeabilityEM() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEM() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEM() const;
            
            virtual scai::lama::Vector<ValueType> const &getVelocityEM();
            virtual scai::lama::Vector<ValueType> const &getVelocityEM() const;
            
            virtual scai::lama::Vector<ValueType> const &getConductivityEMoptical();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMoptical() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMoptical();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMoptical() const;
            virtual ValueType const &getTauDisplacementEM();
            virtual ValueType const &getTauDisplacementEM() const;
            virtual scai::lama::Vector<ValueType> const &getTauConductivityEM() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEM() const;

            virtual scai::lama::Vector<ValueType> const &getPorosity() const;
            virtual scai::lama::Vector<ValueType> const &getSaturation() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEMWater() const;
            virtual scai::lama::Vector<ValueType> const &getRelativeDieletricPeimittivityRockMatrix() const;

            virtual scai::IndexType getNumRelaxationMechanisms() const;
            virtual ValueType getRelaxationFrequency() const;

            virtual void setMagneticPermeabilityEM(scai::lama::Vector<ValueType> const &setMagneticPermeabilityEM);
            virtual void setConductivityEM(scai::lama::Vector<ValueType> const &setConductivityEM);
            virtual void setDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setDielectricPermittivityEM);
            
            virtual void setConductivityEMoptical(scai::lama::Vector<ValueType> const &setConductivityEMoptical);
            virtual void setDielectricPermittivityEMoptical(scai::lama::Vector<ValueType> const &setDielectricPermittivityEMoptical);
            virtual void setTauConductivityEM(scai::lama::Vector<ValueType> const &setTauConductivityEM);
            virtual void setTauDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setTauDielectricPermittivityEM);
            
            virtual void setPorosity(scai::lama::Vector<ValueType> const &setPorosity);
            virtual void setSaturation(scai::lama::Vector<ValueType> const &setSaturation);
            virtual void setConductivityEMWater(scai::lama::Vector<ValueType> const &setConductivityEMWater);
            virtual void setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix); 
            
            ValueType const getDielectricPermittivityVacuum() const; 
            ValueType const getConductivityReference() const; 
            ValueType const getTauDielectricPermittivityReference() const; 
            ValueType const getTauConductivityReference() const; 
            ValueType const getRelativeDielectricPermittivityWater() const; 
            ValueType const getRelativeDielectricPermittivityVacuum() const; 
            
            void setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const);  
            ValueType const getArchie_a() const; 
            ValueType const getArchie_m() const;  
            ValueType const getArchie_n() const;            
            void calcRockMatrixParameter(Configuration::Configuration const &config);
            void calcConductivityReference(ValueType const CenterFrequencyCPML);
            void calcWaveModulusFromPetrophysics();
            void calcPetrophysicsFromWaveModulus();
            
            virtual void setNumRelaxationMechanisms(scai::IndexType const setNumRelaxationMechanisms);
            virtual void setRelaxationFrequency(ValueType const setRelaxationFrequency);
            
            scai::IndexType getParameterisation() const;
            void setParameterisation(scai::IndexType const setParameterisation);
            scai::IndexType getInversionType() const;
            void setInversionType(scai::IndexType const setInversionType);

            /*! \brief Prepare the model parameters for modelling */
            virtual void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) = 0;
            void prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

            virtual void applyThresholds(Configuration::Configuration const &config) = 0;

            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageYZ();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageYZ() const;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageXZ();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageXZ() const;
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageXY();
            virtual scai::lama::Vector<ValueType> const &getInverseMagneticPermeabilityEMAverageXY() const;
            
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageX();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageY();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageZ();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageX();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageY();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageZ();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageX();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageY();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageZ();
            virtual scai::lama::Vector<ValueType> const &getConductivityEMopticalAverageZ() const;            
            
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageX();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageY();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageZ();
            virtual scai::lama::Vector<ValueType> const &getDielectricPermittivityEMopticalAverageZ() const;
            
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageX();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageX() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageY();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageY() const;
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageZ();
            virtual scai::lama::Vector<ValueType> const &getTauDielectricPermittivityEMAverageZ() const;
            
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
            bool dirtyFlagConductivityEMoptical; //!< ==true if conductivityEMoptical has to be recalulated;
            bool dirtyFlagDielectricPermittivityEMoptical; //!< ==true if dielectricPermittivityEMoptical has to be recalulated;

            scai::IndexType parameterisation;
            scai::IndexType inversionType;
            scai::IndexType fileFormat;      //!< 1=mtx 2=lmf

            ValueType aArchie;
            ValueType mArchie;
            ValueType nArchie;   
            ValueType const MagneticPermeabilityVacuum = 1.2566370614e-6;       
            ValueType const DielectricPermittivityVacuum = 8.8541878176e-12;   
            ValueType ConductivityReference;   
            ValueType const TauDielectricPermittivityReference = 0.01; 
            ValueType const RelativeDielectricPermittivityWater = 81; 
            ValueType const RelativeDielectricPermittivityVacuum = 1;
            
            std::string equationType;

            scai::lama::DenseVector<ValueType> velocivityEM;   //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> magneticPermeabilityEM;        //!< Vector storing MagneticPermeabilityEM.
            scai::lama::DenseVector<ValueType> conductivityEM; //!< Vector storing EM-wave velocity.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEM; //!< Vector storing dielectricPermittivityEM.
            
            scai::lama::DenseVector<ValueType> conductivityEMoptical; //!< Vector storing dielectricPermittivityEMoptical for visco-emem modelling.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMoptical; //!< Vector storing dielectricPermittivityEMoptical for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauConductivityEM; //!< Vector storing tauDielectricPermittivityEM for visco-emem modelling.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEM; //!< Vector storing tauDielectricPermittivityEM for visco-emem modelling.
            
            scai::lama::DenseVector<ValueType> porosity; //!< Vector storing porosity.
            scai::lama::DenseVector<ValueType> saturation; //!< Vector storing saturation.
            scai::lama::DenseVector<ValueType> conductivityEMWater; //!< Vector conductivityEMWater.
            scai::lama::DenseVector<ValueType> relativeDieletricPeimittivityRockMatrix; //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityEMAverageYZ; //!< Vector storing inverse averaged magneticPermeabilityEM in yz-plan.
            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityEMAverageXZ; //!< Vector storing inverse averaged magneticPermeabilityEM in xz-plan.
            scai::lama::DenseVector<ValueType> inverseMagneticPermeabilityEMAverageXY; //!< Vector storing inverse averaged magneticPermeabilityEM in xy-plan.
            scai::lama::DenseVector<ValueType> conductivityEMAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> conductivityEMAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> conductivityEMAverageZ; //!< Vector storing averaged modulus in z-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMAverageZ; //!< Vector storing averaged modulus in z-direction.
            
            scai::lama::DenseVector<ValueType> conductivityEMopticalAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> conductivityEMopticalAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> conductivityEMopticalAverageZ; //!< Vector storing averaged modulus in z-direction.          
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMopticalAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMopticalAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> dielectricPermittivityEMopticalAverageZ; //!< Vector storing averaged modulus in z-direction.              
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEMAverageX; //!< Vector storing averaged modulus in x-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEMAverageY; //!< Vector storing averaged modulus in y-direction.
            scai::lama::DenseVector<ValueType> tauDielectricPermittivityEMAverageZ; //!< Vector storing averaged modulus in z-direction.
         
            scai::IndexType numRelaxationMechanisms; //!< Number of relaxation mechanisms
            ValueType relaxationFrequency;           //!< Relaxation Frequency
            ValueType tauDisplacementEM;     // relaxation time of electric displacement

            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value);
            void initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            void calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecMagneticPermeabilityEM, scai::lama::Vector<ValueType> &vecDielectricPermittivityEM);

            void calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecDielectricPermittivityEM, scai::lama::Vector<ValueType> &vecMagneticPermeabilityEM, scai::lama::Vector<ValueType> &vecVelocity);

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

            SparseFormat averageMatrixX;                                 //!< Averaging magneticPermeabilityEM matrix in x-direction
            SparseFormat averageMatrixY;                                 //!< Averaging magneticPermeabilityEM matrix in y-direction
            SparseFormat averageMatrixZ;                                 //!< Averaging magneticPermeabilityEM matrix in z-direction
            SparseFormat averageMatrixXY;                                //!< Average EM-wave Modulus in xy-plane
            SparseFormat averageMatrixXZ;                                //!< Average EM-wave Modulus in xz-plane
            SparseFormat averageMatrixYZ;                                //!< Average EM-wave Modulus in yz-plane

            void calculateInverseAveragedEMparameter(scai::lama::Vector<ValueType> const &vecEMparameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedEMparameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
            void calculateAveragedEMparameter(scai::lama::Vector<ValueType> const &vecEMparameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix);
                         
          private:
            void allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
        };
    }
}
