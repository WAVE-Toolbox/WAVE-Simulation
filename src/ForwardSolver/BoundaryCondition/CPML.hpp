#pragma once

namespace KITGPI {
    
    namespace ForwardSolver {
        
        namespace BoundaryCondition {
            
            //! \brief Abstract class for the calculation of the Absorbing Boundary condition
            template<typename ValueType>
            class CPML
            {
            public:
                
                //! \brief Default constructor
                CPML(){};
                
                //! \brief Default destructor
                ~CPML(){};
                
                
                //! \brief Initializsation of the Absorbing-Coefficient matrix
                /*!
                 *
                 \param dist Distribution of the wavefield
                 \param ctx Context
                 \param NX Total number of grid points in X
                 \param NY Total number of grid points in Y
                 \param NZ Total number of grid points in Z
                 \param BoundaryWidth Width of damping boundary
                 \param useFreeSurface Bool if free surface is in use
                 */
                virtual void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ,ValueType DT, IndexType DH, IndexType BoundaryWidth,bool useFreeSurface,Configuration::PMLVariables<ValueType> &PMLVar)=0;
                                
            };
        } /* end namespace BoundaryCondition  */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */

