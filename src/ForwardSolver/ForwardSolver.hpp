
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

            /*! \brief Running the foward solver
             *
             * Start the forward solver as defined by the given parameters
             *
             \param receiver Configuration of the receivers
             \param sources Configuration of the sources
             \param model Configuration of the modelparameter
             \param wavefield Wavefields for the modelling
             \param derivatives Derivations matrices to calculate the spatial derivatives
             \param NT Total number of time steps
             \param DT Temporal Sampling intervall in seconds
             */
            virtual void run(Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives) = 0;

            /*! \brief Reset PML memory variables if PML is used 
            *
            */
            virtual void resetCPML() = 0;

            /*! \brief Sets Scaling facrors for horizontal updates at the Free Surface
            *
            *
            \param model model parameter object
            */
            virtual void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) = 0;

            ValueType estimateBoundaryMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML);

            virtual ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            /*! \brief Initialitation of the boundary conditions
             *
             *
             \param config Configuration
             \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
             \param derivatives Derivatives matrices
             \param dist Distribution of the wave fields
             \param ctx Context
             */
            virtual void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) = 0;

            void prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType> &FreeSurface, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML);

            /*! \brief Initialitation of the ForwardSolver
             *
             *
             \param config Configuration
             \param derivatives Derivatives matrices
             \param wavefield Wavefields for the modelling
             \param ctx Context
             */
            virtual void initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) = 0;

          protected:
            scai::IndexType useFreeSurface; //!< Indicator which free surface is in use
            bool useDampingBoundary;        //!< Bool if damping boundary is in use
            bool useConvPML;                //!< Bool if CPML is in use
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
