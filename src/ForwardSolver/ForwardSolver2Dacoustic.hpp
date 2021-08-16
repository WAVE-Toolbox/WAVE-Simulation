
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolverSeismic.hpp"

#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2DAcoustic.hpp"
#include "BoundaryCondition/FreeSurface2Dacoustic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 2-D Acoustic forward solver
        template <typename ValueType>
        class FD2Dacoustic : public ForwardSolverSeismic<ValueType>
        {

          public:
            //! Default constructor
            FD2Dacoustic(){};

            //! Default destructor
            ~FD2Dacoustic(){};

            ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/) override;

            void run(Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const & /*model*/, ValueType /*DT*/) override{/*Nothing todo in acoustic modelling*/};

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Dacoustic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML2DAcoustic<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;

            /* Auxiliary Vectors */
            using ForwardSolver<ValueType>::update;
            using ForwardSolver<ValueType>::update_temp;
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
