
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

        //! Class for Modelparameter for visco-emem simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the visco-emem finite-difference simulation.
         */
        template <typename ValueType>
        class ViscoEMEM : public ModelparameterEM<ValueType>
        {
          public:
            //! Default constructor.
            ViscoEMEM() { equationType = "viscoemem"; };

            //! Destructor, releases all allocated resources.
            ~ViscoEMEM(){};

            explicit ViscoEMEM(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType tauConductivityEM_const, ValueType tauDielectricPermittivityEM_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            explicit ViscoEMEM(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            ViscoEMEM(const ViscoEMEM &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType magneticPermeabilityEM_const, ValueType conductivityEM_const, ValueType dielectricPermittivityEM_const, ValueType tauConductivityEM_const, ValueType tauDielectricPermittivityEM_const, scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;

            void initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, ValueType relaxationFrequency_in);

            void write(std::string filename, scai::IndexType fileFormat) const override;

            std::string getEquationType() const;
            
            scai::lama::Vector<ValueType> const &getVelocityEM() override;
            scai::lama::Vector<ValueType> const &getVelocityEM() const override;
            
            /* Getter methods for not requiered parameters */
            
            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void getModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;
            void setModelPerShot(KITGPI::Modelparameter::ModelparameterEM<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth) override;

            void minusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::ViscoEMEM<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> operator+(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> &operator+=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> operator-(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> &operator-=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoEMEM<ValueType> &operator=(KITGPI::Modelparameter::ViscoEMEM<ValueType> const &rhs);

          private:
            void init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates){COMMON_THROWEXCEPTION("variable grid is not implemented in the viscoemem case")};

            void calculateAveraging() override;

            using ModelparameterEM<ValueType>::MagneticPermeabilityVacuum;        // magnetic permeability of free space
            using ModelparameterEM<ValueType>::DielectricPermittivityVacuum;     // dielectric permittivity of free space 
            
            using ModelparameterEM<ValueType>::equationType;

            using ModelparameterEM<ValueType>::dirtyFlagAveraging;
            using ModelparameterEM<ValueType>::dirtyFlagVelocivityEM;            
            using ModelparameterEM<ValueType>::dirtyFlagConductivityEMoptical;
            using ModelparameterEM<ValueType>::dirtyFlagDielectricPermittivityEMoptical;
            
            using ModelparameterEM<ValueType>::velocivityEM;
            using ModelparameterEM<ValueType>::magneticPermeabilityEM;
            using ModelparameterEM<ValueType>::conductivityEM;
            using ModelparameterEM<ValueType>::dielectricPermittivityEM;
            using ModelparameterEM<ValueType>::porosity; //!< Vector storing porosity.
            using ModelparameterEM<ValueType>::saturation; //!< Vector storing saturation.
            using ModelparameterEM<ValueType>::conductivityEMWater;        //!< Vector storing conductivityEMWater.
            using ModelparameterEM<ValueType>::relativeDieletricPeimittivityRockMatrix;   //!< Vector storing relativeDieletricPeimittivityRockMatrix.

            using ModelparameterEM<ValueType>::conductivityEMoptical;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMoptical;
            using ModelparameterEM<ValueType>::tauConductivityEM;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityEM;
            using ModelparameterEM<ValueType>::relaxationFrequency;
            using ModelparameterEM<ValueType>::numRelaxationMechanisms;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;

            void purgeMatrices() override;

            using ModelparameterEM<ValueType>::averageMatrixX;
            using ModelparameterEM<ValueType>::averageMatrixY;
            using ModelparameterEM<ValueType>::averageMatrixZ;
            using ModelparameterEM<ValueType>::averageMatrixXY;
            using ModelparameterEM<ValueType>::averageMatrixXZ;
            using ModelparameterEM<ValueType>::averageMatrixYZ;

            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityEMAverageYZ;
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityEMAverageXZ;
            using ModelparameterEM<ValueType>::inverseMagneticPermeabilityEMAverageXY;
            using ModelparameterEM<ValueType>::conductivityEMopticalAverageY;
            using ModelparameterEM<ValueType>::conductivityEMopticalAverageX;
            using ModelparameterEM<ValueType>::conductivityEMopticalAverageZ;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMopticalAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMopticalAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMopticalAverageZ;
            
            using ModelparameterEM<ValueType>::dielectricPermittivityEMAverageX;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMAverageY;
            using ModelparameterEM<ValueType>::dielectricPermittivityEMAverageZ;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityEMAverageX;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityEMAverageY;
            using ModelparameterEM<ValueType>::tauDielectricPermittivityEMAverageZ;
            using ModelparameterEM<ValueType>::tauDisplacementEM;
            
            /* Not requiered parameters */
        };
    }
}
