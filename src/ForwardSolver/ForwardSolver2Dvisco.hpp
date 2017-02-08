
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Seismogram.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2D.hpp"
#include "BoundaryCondition/FreeSurface2Dvisco.hpp"
#include "Derivatives/Derivatives.hpp"
#include "SourceReceiverImpl/FDTD2Delastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 3-D visco forward solver
        template <typename ValueType>
        class FD2Dvisco : public ForwardSolver<ValueType>
        {

          public:
            //! Default constructor
            FD2Dvisco(){};

            //! Default destructor
            ~FD2Dvisco(){};

            void run(Acquisition::Receivers<ValueType> &receiver, Acquisition::Sources<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType NT, ValueType DT) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Dvisco<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;

            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;

            BoundaryCondition::CPML2D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
