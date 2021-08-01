
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

#include "Modelparameter.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! Class for Modelparameter for viscotmem simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the viscotmem finite-difference simulation.
         */
        template <typename ValueType>
        class ViscoTMEM : public ModelparameterEM<ValueType>
        {
          public:
            //! Default constructor.
            ViscoTMEM() { equationType = "viscotmem"; };

            //! Destructor, releases all allocated resources.
            ~ViscoTMEM(){};

            explicit ViscoTMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit ViscoTMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            explicit ViscoTMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            ViscoTMEM(const ViscoTMEM &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const, ValueType tauElectricConductivity_const, ValueType tauDielectricPermittivity_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            
            void write(std::string filename, scai::IndexType fileFormat) const override;

            std::string getEquationType() const;
            
            scai::lama::Vector<ValueType> const &getVelocityEM() override;
            scai::lama::Vector<ValueType> const &getVelocityEM() const override;
                        
            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;
            void setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) override;

            void minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::ViscoTMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> operator+(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> &operator+=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> operator-(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> &operator-=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoTMEM<ValueType> &operator=(KITGPI::Modelparameter::ViscoTMEM<ValueType> const &rhs);

          private:
            void init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates);
            void calculateAveraging() override;

            using ModelparameterEM<ValueType>::MagneticPermeabilityVacuum;        // magnetic permeability of free space
            using ModelparameterEM<ValueType>::DielectricPermittivityVacuum;     // dielectric permittivity of free space 
            
            using ModelparameterEM<ValueType>::equationType;

            using ModelparameterEM<ValueType>::dirtyFlagAveraging;
            using ModelparameterEM<ValueType>::dirtyFlagVelocivityEM;
            using ModelparameterEM<ValueType>::velocivityEM;
            using ModelparameterEM<ValueType>::magneticPermeability;
            using ModelparameterEM<ValueType>::electricConductivity;
            using ModelparameterEM<ValueType>::dielectricPermittivity;
            using ModelparameterEM<ValueType>::porosity; //!< Vector storing porosity.
            using ModelparameterEM<ValueType>::saturation; //!< Vector storing saturation.
            using ModelparameterEM<ValueType>::electricConductivityWater;        //!< Vector storing electricConductivityWater.
            using ModelparameterEM<ValueType>::relativeDieletricPeimittivityRockMatrix;   //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            using ModelparameterEM<ValueType>::dirtyFlagElectricConductivityOptical;
            using ModelparameterEM<ValueType>::dirtyFlagDielectricPermittivityOptical;
            using ModelparameterEM<ValueType>::electricConductivityOptical;
            using ModelparameterEM<ValueType>::dielectricPermittivityOptical;
            using ModelparameterEM<ValueType>::tauElectricConductivity;
            using ModelparameterEM<ValueType>::tauDielectricPermittivity;
            using ModelparameterEM<ValueType>::relaxationFrequency;
            using ModelparameterEM<ValueType>::numRelaxationMechanisms;
            
            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
            void purgeMatrices() override;

            using ModelparameterEM<ValueType>::averageMatrixX;
            using ModelparameterEM<ValueType>::averageMatrixY;

            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageYZ;
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageXZ;
            
            using ModelparameterEM<ValueType>::tauElectricDisplacement;

            /* Not requiered parameters */
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageXY;
            using ModelparameterEM<ValueType>::electricConductivityAverageX;
            using ModelparameterEM<ValueType>::electricConductivityAverageY;
            using ModelparameterEM<ValueType>::electricConductivityAverageZ;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageZ;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageY;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageX;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageZ;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageZ;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageX;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageY;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageZ;
        };
    }
}
