
#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include <iostream>

#include "ForwardSolver.hpp"

#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2D.hpp"
#include "BoundaryCondition/FreeSurface2Dvisco.hpp"
#include "SourceReceiverImpl/FDTD2Delastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief 2-D visco forward solver
        template <typename ValueType>
        class FD2Dvisco : public ForwardSolver<ValueType>
        {

          public:
            //! Default constructor
            FD2Dvisco(){};

            //! Default destructor
            ~FD2Dvisco(){};

            void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, const Acquisition::AcquisitionGeometry<ValueType> &sources, const Modelparameter::Modelparameter<ValueType> &model, Wavefields::Wavefields<ValueType> &wavefield, const Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType t) override;

            void resetCPML() override;

            void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) override;

            void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates , Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) override;

            void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, const Modelparameter::Modelparameter<ValueType> &model,Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) override;

          private:
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Dvisco<ValueType> FreeSurface; //!< Free Surface boundary condition class
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
            scai::lama::DenseVector<ValueType> update2;
            scai::lama::DenseVector<ValueType> onePlusLtauP;
            scai::lama::DenseVector<ValueType> onePlusLtauS;

            IndexType numRelaxationMechanisms; // = Number of relaxation mechanisms
            ValueType relaxationTime;          // = 1 / ( 2 * Pi * f_relax )
            ValueType inverseRelaxationTime;   // = 1 / relaxationTime
            ValueType viscoCoeff1;             // = 1 - DT / ( 2 * tau_Sigma_l )
            ValueType viscoCoeff2;             // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
            ValueType DThalf;                  // = DT / 2.0
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
