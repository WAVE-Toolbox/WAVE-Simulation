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
                
                //! \brief Getter method for derivative matrix A
                virtual lama::CSRSparseMatrix<ValueType>& getA();
                //! \brief Getter method for derivative matrix B
                virtual lama::CSRSparseMatrix<ValueType>& getB();
                //! \brief Getter method for derivative matrix C
                virtual lama::CSRSparseMatrix<ValueType>& getC();
                //! \brief Getter method for derivative matrix D
                virtual lama::CSRSparseMatrix<ValueType>& getD();
                //! \brief Getter method for derivative matrix E
                virtual lama::CSRSparseMatrix<ValueType>& getE();
                //! \brief Getter method for derivative matrix F
                virtual lama::CSRSparseMatrix<ValueType>& getF();
                
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
                
                lama::CSRSparseMatrix<ValueType> A; //!< Derivative matrix A
                lama::CSRSparseMatrix<ValueType> B; //!< Derivative matrix B
                lama::CSRSparseMatrix<ValueType> C; //!< Derivative matrix C
                lama::CSRSparseMatrix<ValueType> D; //!< Derivative matrix D
                lama::CSRSparseMatrix<ValueType> E; //!< Derivative matrix E
                lama::CSRSparseMatrix<ValueType> F; //!< Derivative matrix F
                
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

//! \brief Getter method for derivative matrix A
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getA(){
    return(A);
}

//! \brief Getter method for derivative matrix B
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getB(){
    return(B);
}

//! \brief Getter method for derivative matrix C
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getC(){
    return(C);
}

//! \brief Getter method for derivative matrix D
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getD(){
    return(D);
}

//! \brief Getter method for derivative matrix E
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getE(){
    return(E);
}

//! \brief Getter method for derivative matrix F
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getF(){
    return(F);
}
