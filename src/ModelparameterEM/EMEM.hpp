
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

        //! Class for Modelparameter for emem simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the emem finite-difference simulation.
         */
        template <typename ValueType>
        class EMEM : public ModelparameterEM<ValueType>
        {
          public:
            //! Default constructor.
            EMEM() { equationType = "emem"; };

            //! Destructor, releases all allocated resources.
            ~EMEM(){};

            explicit EMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const);
            explicit EMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            EMEM(const EMEM &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeability_const, ValueType electricConductivity_const, ValueType dielectricPermittivity_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void write(std::string filename, scai::IndexType fileFormat) const override;

            std::string getEquationType() const;

            /* Getter methods for not required parameters */
            scai::lama::Vector<ValueType> const &getTauElectricConductivity() const override;
            scai::lama::Vector<ValueType> const &getTauDielectricPermittivity() const override;
            
            scai::IndexType getNumRelaxationMechanisms() const override;
            std::vector<ValueType> getRelaxationFrequency() const override;
            
            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;

            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::EMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::EMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::EMEM<ValueType> operator+(KITGPI::Modelparameter::EMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::EMEM<ValueType> &operator+=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::EMEM<ValueType> operator-(KITGPI::Modelparameter::EMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::EMEM<ValueType> &operator-=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::EMEM<ValueType> &operator=(KITGPI::Modelparameter::EMEM<ValueType> const &rhs);

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
            using Modelparameter<ValueType>::reflectivity;
            using Modelparameter<ValueType>::electricConductivityWater;        //!< Vector storing electricConductivityWater.
            using Modelparameter<ValueType>::relativeDieletricPeimittivityRockMatrix;   //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
            void purgeMatrices() override;

            using Modelparameter<ValueType>::averageMatrixX;
            using Modelparameter<ValueType>::averageMatrixY;
            using Modelparameter<ValueType>::averageMatrixZ;
            using Modelparameter<ValueType>::averageMatrixXY;
            using Modelparameter<ValueType>::averageMatrixXZ;
            using Modelparameter<ValueType>::averageMatrixYZ;

            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageYZ;
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXZ;
            using Modelparameter<ValueType>::inverseMagneticPermeabilityAverageXY;
            using Modelparameter<ValueType>::electricConductivityAverageY;
            using Modelparameter<ValueType>::electricConductivityAverageX;
            using Modelparameter<ValueType>::electricConductivityAverageZ;
            using Modelparameter<ValueType>::dielectricPermittivityAverageZ;
            using Modelparameter<ValueType>::dielectricPermittivityAverageY;
            using Modelparameter<ValueType>::dielectricPermittivityAverageX;
            
            /* Not required parameters */
            using Modelparameter<ValueType>::tauElectricConductivity;
            using Modelparameter<ValueType>::tauDielectricPermittivity;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::centerFrequencyCPML;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::tauElectricConductivityAverageX;
            using Modelparameter<ValueType>::tauElectricConductivityAverageY;
            using Modelparameter<ValueType>::tauElectricConductivityAverageZ;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageX;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageY;
            using Modelparameter<ValueType>::tauDielectricPermittivityAverageZ;
            using Modelparameter<ValueType>::relaxationTime;
        };
    }
}
