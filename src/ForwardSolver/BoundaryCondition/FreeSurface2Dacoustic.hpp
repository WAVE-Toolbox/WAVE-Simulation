#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurfaceAcoustic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D acoustic free surface
            template<typename ValueType>
            class FreeSurface2Dacoustic : public FreeSurfaceAcoustic<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface2Dacoustic(){};
                
                //! Default destructor
                ~FreeSurface2Dacoustic(){};
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */

