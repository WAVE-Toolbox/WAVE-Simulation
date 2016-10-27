#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurface.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief Abstract free surface class
            template<typename ValueType>
            class FreeSurface
            {
            public:
                
                //! Default constructor
                FreeSurface(){};
                
                //! Default destructor
                ~FreeSurface(){};
                
                /*! \brief Initialitation of the free surface
                 *
                 *
                 \param dist Distribution of wavefields
                 \param derivatives Derivative class
                 \param NX Number of grid points in X-Direction
                 \param NY Number of grid points in Y-Direction (Depth)
                 \param NZ Number of grid points in Z-Direction
                 */
                virtual void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ)=0;
                
                /*! \brief Getter method for active bool
                 *
                 *
                 */
                virtual bool getActive()=0;
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
