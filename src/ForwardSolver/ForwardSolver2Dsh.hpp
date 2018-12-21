
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML3D.hpp"
#include "BoundaryCondition/FreeSurface2Delastic.hpp"
#include "SourceReceiverImpl/FDTD2Dsh.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 3-D elastic forward solver
        template <typename ValueType>
        class FD2Dsh : public ForwardSolver<ValueType>
        {

          public:
            //! Default constructor
            FD2Dsh(){};

            //! Default destructor
            ~FD2Dsh(){};

            void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, const Acquisition::AcquisitionGeometry<ValueType> &sources, const Modelparameter::Modelparameter<ValueType> &model, Wavefields::Wavefields<ValueType> &wavefield, const Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType t) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates const &modelCoordinates , Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model,Acquisition::Coordinates const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/) override;
	    
          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML3D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;
	    
	                /* Auxiliary Vectors */
            scai::lama::DenseVector<ValueType> update;
            scai::lama::DenseVector<ValueType> update_temp;

        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
