
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2D.hpp"
#include "BoundaryCondition/FreeSurface2Demem.hpp"
#include "SourceReceiverImpl/FDTD2Demem.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 3-D emem forward solver
        template <typename ValueType>
        class FD2Demem : public ForwardSolverEM<ValueType>
        {

          public:
            //! Default constructor
            FD2Demem(){};

            //! Default destructor
            ~FD2Demem(){};

            ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void run(Acquisition::AcquisitionGeometryEM<ValueType> &receiver, Acquisition::AcquisitionGeometryEM<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) override;
            void runAdjoint(Acquisition::AcquisitionGeometryEM<ValueType> &receiver, Acquisition::AcquisitionGeometryEM<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::ModelparameterEM<ValueType> const &model, ValueType /*DT*/) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Wavefields::WavefieldsEM<ValueType> &wavefield, Modelparameter::ModelparameterEM<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Demem<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolverEM<ValueType>::useFreeSurface;

            BoundaryCondition::ABSEM2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolverEM<ValueType>::useDampingBoundary;

            BoundaryCondition::CPMLEM2D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolverEM<ValueType>::useConvPML;
            
            /* Auxiliary Vectors */
            using ForwardSolverEM<ValueType>::update;
            using ForwardSolverEM<ValueType>::update_temp;
            using ForwardSolverEM<ValueType>::DT_temp;
            using ForwardSolverEM<ValueType>::CaAverageX;
            using ForwardSolverEM<ValueType>::CaAverageY;
            using ForwardSolverEM<ValueType>::CbAverageX;
            using ForwardSolverEM<ValueType>::CbAverageY;
            
            /* Not requiered coefficients */
            using ForwardSolverEM<ValueType>::CaAverageZ;
            using ForwardSolverEM<ValueType>::CbAverageZ;
            using ForwardSolverEM<ValueType>::Cc;
            using ForwardSolverEM<ValueType>::CdAverageX;
            using ForwardSolverEM<ValueType>::CdAverageY;
            using ForwardSolverEM<ValueType>::CdAverageZ;
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
