#pragma once

#include "../ForwardSolver/ForwardSolver.hpp"
#include "BoundaryCondition/FreeSurfaceEM.hpp"

namespace KITGPI
{

    //! \brief ForwardSolver namespace
    namespace ForwardSolver
    {

        //! \brief Abstract class for forward solver
        template <typename ValueType>
        class ForwardSolverEM : public ForwardSolver<ValueType>
        {
          public:
            //! Default constructor
            ForwardSolverEM(){};

            //! Default destructor
            ~ForwardSolverEM(){};

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
            scai::lama::DenseVector<ValueType> getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) override;
            scai::lama::DenseVector<ValueType> getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) override;
            scai::lama::DenseVector<ValueType> getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT);
            ValueType getCc(ValueType tauElectricDisplacement, ValueType DT) override;
            scai::lama::DenseVector<ValueType> getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, ValueType tauElectricDisplacement, ValueType DT) override;
                  
            ValueType DT_temp; //!< Vector storing DT.
            scai::lama::DenseVector<ValueType> CaAverageX; //!< Vector storing averaged coefficient Ca in x-direction.
            scai::lama::DenseVector<ValueType> CaAverageY; //!< Vector storing averaged coefficient Ca in y-direction.
            scai::lama::DenseVector<ValueType> CaAverageZ; //!< Vector storing averaged coefficient Ca in z-direction.
            scai::lama::DenseVector<ValueType> CbAverageX; //!< Vector storing averaged coefficient Cb in x-direction.
            scai::lama::DenseVector<ValueType> CbAverageY; //!< Vector storing averaged coefficient Cb in y-direction.
            scai::lama::DenseVector<ValueType> CbAverageZ; //!< Vector storing averaged coefficient Cb in z-direction.            
            ValueType Cc; //!< Vector storing averaged coefficient Cc.
            scai::lama::DenseVector<ValueType> CdAverageX; //!< Vector storing averaged coefficient Cd in x-direction.
            scai::lama::DenseVector<ValueType> CdAverageY; //!< Vector storing averaged coefficient Cd in y-direction.
            scai::lama::DenseVector<ValueType> CdAverageZ; //!< Vector storing averaged coefficient Cd in z-direction.
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
