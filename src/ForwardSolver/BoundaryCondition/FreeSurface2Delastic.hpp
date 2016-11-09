#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurfaceElastic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D elastic free surface
            template<typename ValueType>
            class FreeSurface2Delastic : public FreeSurfaceElastic<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface2Delastic(){};
                
                //! Default destructor
                ~FreeSurface2Delastic(){};
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy);

            private:
                
                using FreeSurfaceElastic<ValueType>::setSurfaceZero;
                using FreeSurfaceElastic<ValueType>::scaleHorizontalUpdate;
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Apply free surface condition during time stepping for 2D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy){
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    
    Sxx +=sumHorizonalDerivative;
    
    Syy.scale(setSurfaceZero);
    
}

