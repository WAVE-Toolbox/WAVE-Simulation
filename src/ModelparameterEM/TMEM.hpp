
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

        //! Class for Modelparameter for tmem simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the tmem finite-difference simulation.
         */
        template <typename ValueType>
        class TMEM : public ModelparameterEM<ValueType>
        {
          public:
            //! Default constructor.
            TMEM() { equationType = "tmem"; };

            //! Destructor, releases all allocated resources.
            ~TMEM(){};

            explicit TMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit TMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const);
            explicit TMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            TMEM(const TMEM &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void write(std::string filename, scai::IndexType fileFormat) const override;

            std::string getEquationType() const;
            
            /* Getter methods for not requiered parameters */
            scai::lama::Vector<ValueType> const &getElectricConductivityOptical() override;
            scai::lama::Vector<ValueType> const &getDielectricPermittivityOptical() override;
            ValueType const &getTauElectricDisplacement() override;
            ValueType const &getTauElectricDisplacement() const override;
            scai::lama::Vector<ValueType> const &getTauElectricConductivity() const override;
            scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const override;
            
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;
            
            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;
            void setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) override;

            void minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::TMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::TMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::TMEM<ValueType> operator+(KITGPI::Modelparameter::TMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::TMEM<ValueType> &operator+=(KITGPI::Modelparameter::TMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::TMEM<ValueType> operator-(KITGPI::Modelparameter::TMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::TMEM<ValueType> &operator-=(KITGPI::Modelparameter::TMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::TMEM<ValueType> &operator=(KITGPI::Modelparameter::TMEM<ValueType> const &rhs);

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
            using ModelparameterEM<ValueType>::reflectivity; //!< Vector storing reflectivity.
            using ModelparameterEM<ValueType>::electricConductivityWater;        //!< Vector storing electricConductivityWater.
            using ModelparameterEM<ValueType>::relativeDieletricPeimittivityRockMatrix;   //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
            void purgeMatrices() override;

            using ModelparameterEM<ValueType>::averageMatrixX;
            using ModelparameterEM<ValueType>::averageMatrixY;

            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageYZ;
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageXZ;            

            /* Not requiered parameters */
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityAverageXY;
            using ModelparameterEM<ValueType>::electricConductivityAverageX;
            using ModelparameterEM<ValueType>::electricConductivityAverageY;
            using ModelparameterEM<ValueType>::electricConductivityAverageZ;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityAverageZ;
            using ModelparameterEM<ValueType>::dirtyFlagElectricConductivityOptical;
            using ModelparameterEM<ValueType>::dirtyFlagDielectricPermittivityOptical;
            using ModelparameterEM<ValueType>::electricConductivityOptical;
            using ModelparameterEM<ValueType>::dielectricPermittivityOptical;
            using ModelparameterEM<ValueType>::tauElectricConductivity;
            using ModelparameterEM<ValueType>::tauDielectricPermittivity;
            using ModelparameterEM<ValueType>::relaxationFrequency;
            using ModelparameterEM<ValueType>::numRelaxationMechanisms;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageX;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageY;
            using ModelparameterEM<ValueType>::electricConductivityOpticalAverageZ;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityOpticalAverageZ;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageX;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageY;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityAverageZ;
            using ModelparameterEM<ValueType>::tauElectricDisplacement;
        };
    }
}
