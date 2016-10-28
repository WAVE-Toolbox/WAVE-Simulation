#pragma once

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief Derivatives namespace
        namespace Derivatives {
            
            //! \brief Abstract class for the calculation of the Derivatives matrices
            template<typename ValueType>
            class Derivatives
            {
            public:
                
                //! \brief Default constructor
                Derivatives():spatialFDorder(0){};
                
                //! \brief Default destructor
                ~Derivatives(){};
                
                //! \brief Getter method for derivative matrix Dxf
                virtual lama::CSRSparseMatrix<ValueType>& getDxf();
                //! \brief Getter method for derivative matrix Dyf
                virtual lama::CSRSparseMatrix<ValueType>& getDyf();
                //! \brief Getter method for derivative matrix Dzf
                virtual lama::CSRSparseMatrix<ValueType>& getDzf();
                //! \brief Getter method for derivative matrix Dxb
                virtual lama::CSRSparseMatrix<ValueType>& getDxb();
                //! \brief Getter method for derivative matrix Dyb
                virtual lama::CSRSparseMatrix<ValueType>& getDyb();
                //! \brief Getter method for derivative matrix Dzb
                virtual lama::CSRSparseMatrix<ValueType>& getDzb();
                
                //! \brief Initializsation of the derivative matrices
                /*!
                 *
                 \param dist Distribution of the wavefield
                 \param ctx Context
                 \param NX Total number of grid points in X
                 \param NY Total number of grid points in Y
                 \param NZ Total number of grid points in Z
                 \param DH Grid spacing (equidistant)
                 \param DT Temporal sampling interval
                 \param spatialFDorder FD-order of spatial stencils
                 \param comm Communicator
                 */
                virtual void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, dmemo::CommunicatorPtr comm )=0;
                
                void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
                
                IndexType getSpatialFDorder();
                
            protected:
                
                lama::CSRSparseMatrix<ValueType> Dxf; //!< Derivative matrix Dxf
                lama::CSRSparseMatrix<ValueType> Dyf; //!< Derivative matrix Dyf
                lama::CSRSparseMatrix<ValueType> Dzf; //!< Derivative matrix Dzf
                lama::CSRSparseMatrix<ValueType> Dxb; //!< Derivative matrix Dxb
                lama::CSRSparseMatrix<ValueType> Dyb; //!< Derivative matrix Dyb
                lama::CSRSparseMatrix<ValueType> Dzb; //!< Derivative matrix Dzb
                
                IndexType spatialFDorder; //!< FD-Order of spatial derivative stencils
                
            };
        } /* end namespace Derivatives */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */

//! \brief Getter method for the spatial FD-order
template<typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getSpatialFDorder(){
    return(spatialFDorder);
}

//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    initializeMatrices(dist,ctx,config.getNX(), config.getNY(), config.getNZ(), config.getDH(), config.getDT(), config.getSpatialFDorder(), comm);
}

//! \brief Getter method for derivative matrix Dxf
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxf(){
    return(Dxf);
}

//! \brief Getter method for derivative matrix Dyf
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyf(){
    return(Dyf);
}

//! \brief Getter method for derivative matrix Dzf
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzf(){
    return(Dzf);
}

//! \brief Getter method for derivative matrix Dxb
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxb(){
    return(Dxb);
}

//! \brief Getter method for derivative matrix Dyb
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyb(){
    return(Dyb);
}

//! \brief Getter method for derivative matrix Dzb
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzb(){
    return(Dzb);
}
