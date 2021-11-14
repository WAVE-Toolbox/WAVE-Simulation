
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolverSeismic.hpp"

#include "BoundaryCondition/ABS3D.hpp"
#include "BoundaryCondition/CPML3D.hpp"
#include "BoundaryCondition/FreeSurface3Dviscoelastic.hpp"
#include "SourceReceiverImpl/FDTD3Delastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 3-D viscoelastic forward solver
        template <typename ValueType>
        class FD3Dviscoelastic : public ForwardSolverSeismic<ValueType>
        {

          public:
            //! Default constructor
            FD3Dviscoelastic(){};

            //! Default destructor
            ~FD3Dviscoelastic(){};

            ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Dviscoelastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS3D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML3D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;

            /* Auxiliary Vectors */
            using ForwardSolver<ValueType>::update;
            using ForwardSolver<ValueType>::update_temp;
            scai::lama::DenseVector<ValueType> vxx;
            scai::lama::DenseVector<ValueType> vyy;
            scai::lama::DenseVector<ValueType> vzz;
            scai::lama::DenseVector<ValueType> update2;
            scai::lama::DenseVector<ValueType> onePlusLtauP;
            scai::lama::DenseVector<ValueType> onePlusLtauS;

            IndexType numRelaxationMechanisms; // = Number of relaxation mechanisms
            std::vector<ValueType> relaxationTime;          // = 1 / ( 2 * Pi * f_relax )
            std::vector<ValueType> inverseRelaxationTime;   // = 1 / relaxationTime
            std::vector<ValueType> viscoCoeff1;             // = 1 - DT / ( 2 * tau_Sigma_l )
            std::vector<ValueType> viscoCoeff2;             // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
            ValueType DThalf;                  // = DT / 2.0
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
