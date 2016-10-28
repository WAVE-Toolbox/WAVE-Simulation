#pragma once

namespace KITGPI {
    
    namespace ForwardSolver {
        
        namespace BoundaryCondition {
            
            //! \brief Abstract class for the calculation of the Absorbing Boundary condition
            template<typename ValueType>
            class ABS
            {
            public:
                
                //! \brief Default constructor
                ABS(){};
                
                //! \brief Default destructor
                ~ABS(){};
                
                
                //! \brief Initializsation of the Absorbing-Coefficient matrix
                /*!
                 *
                 \param dist Distribution of the wavefield
                 \param ctx Context
                 \param NX Total number of grid points in X
                 \param NY Total number of grid points in Y
                 \param NZ Total number of grid points in Z
                 \param comm Communicator
                 */
                virtual void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )=0;
                
                virtual void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
                
            };
        } /* end namespace BoundaryCondition  */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS<ValueType>::init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    init(dist,ctx,config.getNX(), config.getNY(), config.getNZ(),comm);
}
