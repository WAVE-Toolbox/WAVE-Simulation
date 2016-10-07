
#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "../Derivatives/Derivatives.hpp"

namespace KITGPI {
    
    //! \brief ForwardSolver namespace
    namespace ForwardSolver {
        
        //! \brief Abstract class for forward solver
        template<typename ValueType>
        class ForwardSolver
        {
        public:
            
            //! Default constructor
            ForwardSolver(){};
            
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
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
