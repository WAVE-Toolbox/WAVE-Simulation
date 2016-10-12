#pragma once

namespace KITGPI {
    
    //! \brief Derivatives namespace
    namespace Derivatives {
        
        //! \brief Abstract class for the calculation of the Derivatives matrices
        template<typename ValueType>
        class Derivatives
        {
        public:
            
            //! \brief Default constructor
            Derivatives(){};

            //! \brief Default destructor
            ~Derivatives(){};
            
            //! \brief Getter method for derivative matrix A
            virtual lama::CSRSparseMatrix<ValueType>& getA()=0;
            //! \brief Getter method for derivative matrix B
            virtual lama::CSRSparseMatrix<ValueType>& getB()=0;
            //! \brief Getter method for derivative matrix C
            virtual lama::CSRSparseMatrix<ValueType>& getC()=0;
            //! \brief Getter method for derivative matrix D
            virtual lama::CSRSparseMatrix<ValueType>& getD()=0;
            //! \brief Getter method for derivative matrix E
            virtual lama::CSRSparseMatrix<ValueType>& getE()=0;
            //! \brief Getter method for derivative matrix F
            virtual lama::CSRSparseMatrix<ValueType>& getF()=0;
            
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
             \param comm Communicator
             */
            virtual void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )=0;
            
            virtual void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );
            
        };
    } /* end namespace Derivatives */
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
void KITGPI::Derivatives::Derivatives<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
    initializeMatrices(dist,ctx,config.getNX(), config.getNY(), config.getNZ(), config.getDH(), config.getDT(),comm);
}
