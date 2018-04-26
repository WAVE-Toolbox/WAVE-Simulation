
#pragma once
#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>

#include "../Acquisition/AcquisitionGeometry.hpp"

#include "../Common/HostPrint.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "Derivatives/Derivatives.hpp"

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
            virtual void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType TStart, scai::IndexType TEnd, ValueType DT) = 0;

            /*! \brief Initialitation of the boundary conditions
             *
             *
             \param config Configuration
             \param derivatives Derivatives matrices
             \param dist Distribution of the wave fields
             \param ctx Context
             */
            virtual void prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) = 0;

          protected:
            scai::IndexType useFreeSurface;     //!< Indicator which free surface is in use
            bool useDampingBoundary; //!< Bool if damping boundary is in use
            bool useConvPML;         //!< Bool if CPML is in use
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
