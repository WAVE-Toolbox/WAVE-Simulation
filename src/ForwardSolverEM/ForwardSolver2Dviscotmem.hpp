
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolverEM.hpp"

#include "../ForwardSolver/BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPMLEM2D.hpp"
#include "SourceReceiverImpl/FDTD2Dtmem.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 2-D tmem forward solver
        template <typename ValueType>
        class FD2Dviscotmem : public ForwardSolverEM<ValueType>
        {

          public:
            //! Default constructor
            FD2Dviscotmem(){};

            //! Default destructor
            ~FD2Dviscotmem(){};

            ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurfaceEM<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPMLEM2D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;

            /* Auxiliary Vectors */
            using ForwardSolver<ValueType>::update;
            using ForwardSolver<ValueType>::update_temp;
            using ForwardSolverEM<ValueType>::DT_temp;
            using ForwardSolverEM<ValueType>::CaAverageZ;
            using ForwardSolverEM<ValueType>::CbAverageZ;
            using ForwardSolverEM<ValueType>::Cc;
            using ForwardSolverEM<ValueType>::CdAverageZ;
            
            /* Not required coefficients */
            using ForwardSolverEM<ValueType>::CaAverageX;
            using ForwardSolverEM<ValueType>::CaAverageY;
            using ForwardSolverEM<ValueType>::CbAverageX;
            using ForwardSolverEM<ValueType>::CbAverageY;
            using ForwardSolverEM<ValueType>::CdAverageX;
            using ForwardSolverEM<ValueType>::CdAverageY;
            
            IndexType numRelaxationMechanisms; // = Number of relaxation mechanisms
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
