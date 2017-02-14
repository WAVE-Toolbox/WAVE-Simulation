
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "BoundaryCondition/ABS3D.hpp"
#include "BoundaryCondition/CPML3DAcoustic.hpp"
#include "BoundaryCondition/FreeSurface3Dacoustic.hpp"
#include "SourceReceiverImpl/FDTD3Dacoustic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 3-D Acoustic forward solver
        template <typename ValueType>
        class FD3Dacoustic : public ForwardSolver<ValueType>
        {

          public:
            //! Default constructor
            FD3Dacoustic(){};

            //! Default destructor
            ~FD3Dacoustic(){};

            void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType NT, ValueType DT) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Dacoustic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS3D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML3DAcoustic<ValueType> ConvPML; //!< CPML boundary condition class
            using ForwardSolver<ValueType>::useConvPML;
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
