
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
            
            //! Default constructor
            ForwardSolver():useFreeSurface(false){};
            
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
             \param comm Communicator
             */
            virtual void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm)=0;
            
            
            /*! \brief Saving seismograms during time stepping
             *
             * THIS METHOD IS CALLED DURING TIME STEPPING
             * DO NOT WASTE RUNTIME HERE
             *
             \param wavefield Wavefields
             \param NT Total number of time steps
             \param t Current time step
             */
            virtual void gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)=0;
            
            /*! \brief Appling the sources to the wavefield
             *
             * THIS METHOD IS CALLED DURING TIME STEPPING
             * DO NOT WASTE RUNTIME HERE
             *
             \param sources Sources to apply
             \param wavefield Wavefields
             \param NT Total number of time steps
             \param t Current time step
             */
            virtual void applySource(Acquisition::Sources<ValueType>& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)=0;
            
            /*! \brief Initialitation of the boundary conditions
             *
             *
             \param config Configuration
             \param derivatives Derivatives matrices
             \param dist Distribution of the wave fields
             */
            virtual void prepareBoundaryConditions(Configuration::Configuration<ValueType> config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist)=0;
            
            Acquisition::Seismogram<ValueType> seismogram; //!< Storage of seismogram data
            
        protected:
            
            bool useFreeSurface; //!< Bool if free surface is in use

            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
