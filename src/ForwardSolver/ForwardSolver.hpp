#pragma once
#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>

#include "../Acquisition/AcquisitionGeometry.hpp"

#include "../Common/HostPrint.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "Derivatives/Derivatives.hpp"
#include "SourceReceiverImpl/SourceReceiverImplFactory.hpp"

#include "BoundaryCondition/ABS.hpp"
#include "BoundaryCondition/CPML.hpp"
#include "BoundaryCondition/FreeSurface.hpp"

namespace KITGPI
{

    //! \brief ForwardSolver namespace
    namespace ForwardSolver
    {

        //! \brief Abstract class for forward solver
        template <typename ValueType>
        class ForwardSolver
        {
          public:
            //! \brief Declare ForwardSolver pointer
            typedef std::shared_ptr<ForwardSolver<ValueType>> ForwardSolverPtr;

            //! Default constructor
            ForwardSolver() : useFreeSurface(false), useDampingBoundary(false), useConvPML(false){};

            //! Default destructor
            ~ForwardSolver(){};

            virtual void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) = 0;

            virtual void resetCPML() = 0;

            virtual void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) = 0;

            ValueType estimateBoundaryMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML);

            virtual ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            virtual void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) = 0;

            void prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType> &FreeSurface, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML);

            virtual void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) = 0;

          protected:
            /* Common */
            scai::IndexType useFreeSurface; //!< Indicator which free surface is in use
            bool useDampingBoundary;        //!< Bool if damping boundary is in use
            bool useConvPML;                //!< Bool if CPML is in use
            
            /* Auxiliary Vectors */
            scai::lama::DenseVector<ValueType> update;
            scai::lama::DenseVector<ValueType> update_temp;
                        
            /* EM */
            virtual scai::lama::DenseVector<ValueType> getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) = 0;
            virtual scai::lama::DenseVector<ValueType> getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) = 0;
            virtual scai::lama::DenseVector<ValueType> getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT) = 0;
            virtual std::vector<ValueType> getCc(scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT) = 0;
            virtual std::vector<scai::lama::DenseVector<ValueType>> getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT) = 0;                  
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
