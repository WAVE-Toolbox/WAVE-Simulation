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
                                
                void setFDCoef(IndexType spFDo);
                
            protected:
                
                lama::CSRSparseMatrix<ValueType> Dxf; //!< Derivative matrix Dxf
                lama::CSRSparseMatrix<ValueType> Dyf; //!< Derivative matrix Dyf
                lama::CSRSparseMatrix<ValueType> Dzf; //!< Derivative matrix Dzf
                lama::CSRSparseMatrix<ValueType> Dxb; //!< Derivative matrix Dxb
                lama::CSRSparseMatrix<ValueType> Dyb; //!< Derivative matrix Dyb
                lama::CSRSparseMatrix<ValueType> Dzb; //!< Derivative matrix Dzb
                
                IndexType spatialFDorder; //!< FD-Order of spatial derivative stencils
                
                scai::hmemo::HArray<ValueType> FDCoef_f;
                scai::hmemo::HArray<ValueType> FDCoef_b;
                
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

//! \brief Set FD coefficients for each order
template<typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDCoef(IndexType spFDo){
    FDCoef_f.resize(spFDo/2);
    FDCoef_b.resize(spFDo/2);
    scai::hmemo::WriteAccess<ValueType> write_FDCoef_f(FDCoef_f);
    scai::hmemo::WriteAccess<ValueType> write_FDCoef_b(FDCoef_b);
    
    switch(spFDo)
    {
        case 2:
            write_FDCoef_f[0]=1.0;
            write_FDCoef_b[0]=-1.0;
            break;
        case 4:
            write_FDCoef_f[1]=-1.0/24.0;
            write_FDCoef_f[0]=9.0/8.0;
            write_FDCoef_b[0]=-9.0/8.0;
            write_FDCoef_b[1]=1.0/24.0;
            break;
        case 6:
            write_FDCoef_f[2]=3.0/640.0;
            write_FDCoef_f[1]=-25.0/384.0;
            write_FDCoef_f[0]=75.0/64.0;
            write_FDCoef_b[0]=-75.0/64.0;
            write_FDCoef_b[1]=25.0/384.0;
            write_FDCoef_b[2]=-3.0/640.0;
            break;
        case 8:
            write_FDCoef_f[3]=-5.0/7168.0;
            write_FDCoef_f[2]=49.0/5120.0;
            write_FDCoef_f[1]=-245.0/3072.0;
            write_FDCoef_f[0]=1225.0/1024.0;
            write_FDCoef_b[0]=-1225.0/1024.0;
            write_FDCoef_b[1]=245.0/3072.0;
            write_FDCoef_b[2]=-49.0/5120.0;
            write_FDCoef_b[3]=5.0/7168.0;
            break;
        case 10:
            write_FDCoef_f[4]=8756999275442633.0/73786976294838206464.0;
            write_FDCoef_f[3]=-8142668969129685.0/4611686018427387904.0;
            write_FDCoef_f[2]=567.0/40960.0;
            write_FDCoef_f[1]=-735.0/8192.0;
            write_FDCoef_f[0]=19845.0/16384.0;
            write_FDCoef_b[0]=-19845.0/16384.0;
            write_FDCoef_b[1]=735.0/8192.0;
            write_FDCoef_b[2]=-567.0/40960.0;
            write_FDCoef_b[3]=8142668969129685.0/4611686018427387904.0;
            write_FDCoef_b[4]=-8756999275442633.0/73786976294838206464.0;
            break;
        case 12:
            write_FDCoef_f[5]=-6448335830095439.0/295147905179352825856.0;
            write_FDCoef_f[4]=1655620175512543.0/4611686018427387904.0;
            write_FDCoef_f[3]=-6842103786556949.0/2305843009213693952.0;
            write_FDCoef_f[2]=628618285389933.0/36028797018963968.0;
            write_FDCoef_f[1]=-436540475965291.0/4503599627370496.0;
            write_FDCoef_f[0]=2750204998582123.0/2251799813685248.0;
            write_FDCoef_b[0]=-2750204998582123.0/2251799813685248.0;
            write_FDCoef_b[1]=436540475965291.0/4503599627370496.0;
            write_FDCoef_b[2]=-628618285389933.0/36028797018963968.0;
            write_FDCoef_b[3]=6842103786556949.0/2305843009213693952.0;
            write_FDCoef_b[4]=-1655620175512543.0/4611686018427387904.0;
            write_FDCoef_b[5]=6448335830095439.0/295147905179352825856.0;
            break;
        default:
            COMMON_THROWEXCEPTION(" Unkown spatialFDorder value.");
            break;
    }
    write_FDCoef_f.release();
    write_FDCoef_b.release();
}
