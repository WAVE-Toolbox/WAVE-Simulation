
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

#include "ModelparameterEM.hpp"

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

            void getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;
            void setModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) override;

            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;

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

            using Modelparameter<ValueType>::MagneticPermeabilityVacuum;        // magnetic permeability of free space
            using Modelparameter<ValueType>::DielectricPermittivityVacuum;     // dielectric permittivity of free space 
            
            using Modelparameter<ValueType>::equationType;

            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::dirtyFlagVelocivityEM;
            using Modelparameter<ValueType>::velocivityEM;
            using Modelparameter<ValueType>::magneticPermeability;
            using Modelparameter<ValueType>::electricConductivity;
            using Modelparameter<ValueType>::dielectricPermittivity;
            using Modelparameter<ValueType>::porosity; //!< Vector storing porosity.
            using Modelparameter<ValueType>::saturation; //!< Vector storing saturation.
            using Modelparameter<ValueType>::reflectivity; //!< Vector storing reflectivity.
            using Modelparameter<ValueType>::electricConductivityWater;        //!< Vector storing electricConductivityWater.
            using Modelparameter<ValueType>::relativeDieletricPeimittivityRockMatrix;   //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
            void purgeMatrices() override;

            using Modelparameter<ValueType>::averageMatrixX;
            using Modelparameter<ValueType>::averageMatrixY;

            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageYZ;
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXZ;            

            /* Not requiered parameters */
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXY;
            using Modelparameter<ValueType>::electricConductivityAverageX;
            using Modelparameter<ValueType>::electricConductivityAverageY;
            using Modelparameter<ValueType>::electricConductivityAverageZ;
            using Modelparameter<ValueType>::dielectricPermittivityAverageX;
            using Modelparameter<ValueType>::dielectricPermittivityAverageY;
            using Modelparameter<ValueType>::dielectricPermittivityAverageZ;
            using Modelparameter<ValueType>::dirtyFlagElectricConductivityOptical;
            using Modelparameter<ValueType>::dirtyFlagDielectricPermittivityOptical;
            using Modelparameter<ValueType>::electricConductivityOptical;
            using Modelparameter<ValueType>::dielectricPermittivityOptical;
            using Modelparameter<ValueType>::tauElectricConductivity;
            using Modelparameter<ValueType>::tauDielectricPermittivity;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::electricConductivityOpticalAverageX;
            using Modelparameter<ValueType>::electricConductivityOpticalAverageY;
            using Modelparameter<ValueType>::electricConductivityOpticalAverageZ;
            using Modelparameter<ValueType>::dielectricPermittivityOpticalAverageX;
            using Modelparameter<ValueType>::dielectricPermittivityOpticalAverageY;
            using Modelparameter<ValueType>::dielectricPermittivityOpticalAverageZ;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageX;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageY;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageZ;
            using Modelparameter<ValueType>::tauElectricDisplacement;
        };
    }
}
