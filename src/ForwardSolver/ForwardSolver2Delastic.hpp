
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2D.hpp"
#include "BoundaryCondition/FreeSurface2Delastic.hpp"
#include "SourceReceiverImpl/FDTD2Delastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 2-D elastic forward solver
        template <typename ValueType>
        class FD2Delastic : public ForwardSolver<ValueType>
        {

          public:
            //! Default constructor
            FD2Delastic(){};

            //! Default destructor
            ~FD2Delastic(){};

            ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void run(Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML2D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;

            /* Auxiliary Vectors */
            scai::lama::DenseVector<ValueType> update;
            scai::lama::DenseVector<ValueType> update_temp;
            scai::lama::DenseVector<ValueType> vxx;
            scai::lama::DenseVector<ValueType> vyy;
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
