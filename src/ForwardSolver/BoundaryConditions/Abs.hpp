#pragma once

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief Abs namespace
        namespace Abs {
            
            //! \brief Abstract class for the calculation of the Absorbing Boundary condition
            template<typename ValueType>
            class Abs
            {
            public:
                
                //! \brief Default constructor
                Abs(){};
                
                //! \brief Default destructor
                ~Abs(){};
                
                //! \brief Getter method for the absorbing-coefficient matrix 
                virtual lama::CSRSparseMatrix<ValueType>& getAbsCoeff()=0;
                
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
                virtual void initializeMatrix(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )=0;
                
                virtual void initializeMatrix(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
                
            };
        } /* end namespace Abs  */
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
void KITGPI::ForwardSolver::Abs::Abs<ValueType>::initializeMatrix(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    initializeMatrix(dist,ctx,config.getNX(), config.getNY(), config.getNZ(),comm);
}
