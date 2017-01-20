
#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "Derivatives/Derivatives.hpp"
#include "../Common/HostPrint.hpp"

namespace KITGPI {
    
    //! \brief ForwardSolver namespace
    namespace ForwardSolver {
        
        //! \brief Abstract class for forward solver
        template<typename ValueType>
        class ForwardSolver
        {
        public:
            
            typedef std::shared_ptr<ForwardSolver<ValueType>> ForwardSolverPtr;
            
            //! Default constructor
            ForwardSolver():useFreeSurface(false),useDampingBoundary(false){};
            
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
            virtual void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType> const& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT)=0;
            
            /*! \brief Initialitation of the boundary conditions
             *
             *
             \param config Configuration
             \param derivatives Derivatives matrices
             \param dist Distribution of the wave fields
             \param ctx Context
             */
            virtual void prepareBoundaryConditions(Configuration::Configuration const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx)=0;
                        
        protected:
            
            bool useFreeSurface; //!< Bool if free surface is in use
            bool useDampingBoundary; //!< Bool if damping boundary is in use

            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
