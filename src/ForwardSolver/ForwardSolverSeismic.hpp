#pragma once

#include "ForwardSolver.hpp"

namespace KITGPI
{

    //! \brief ForwardSolver namespace
    namespace ForwardSolver
    {

        //! \brief Abstract class for forward solver
        template <typename ValueType>
        class ForwardSolverSeismic : public ForwardSolver<ValueType>
        {
          public:
            //! Default constructor
            ForwardSolverSeismic(){};

            //! Default destructor
            ~ForwardSolverSeismic(){};

            virtual void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) = 0;

            virtual void resetCPML() = 0;

            virtual void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) = 0;

            virtual ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) = 0;

            virtual void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) = 0;

          protected:
            /* Common */
            using ForwardSolver<ValueType>::useFreeSurface; //!< Indicator which free surface is in use
            using ForwardSolver<ValueType>::useDampingBoundary;        //!< Bool if damping boundary is in use
            using ForwardSolver<ValueType>::useConvPML;                //!< Bool if CPML is in use       
            using ForwardSolver<ValueType>::update;
            using ForwardSolver<ValueType>::update_temp;      

            /* EM */
            virtual scai::lama::DenseVector<ValueType> getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) override;
            virtual scai::lama::DenseVector<ValueType> getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) override;
            virtual scai::lama::DenseVector<ValueType> getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) override;
            virtual std::vector<ValueType> getCc(scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT) override;
            virtual std::vector<scai::lama::DenseVector<ValueType>> getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT) override;                  
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
