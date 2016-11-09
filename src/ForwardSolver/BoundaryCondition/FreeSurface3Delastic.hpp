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
            class FreeSurface3Delastic : public FreeSurfaceElastic<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface3Delastic(){};
                
                //! Default destructor
                ~FreeSurface3Delastic(){};
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz);

            private:
                
                using FreeSurfaceElastic<ValueType>::setSurfaceZero;
                using FreeSurfaceElastic<ValueType>::scaleHorizontalUpdate;
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Apply free surface condition during time stepping for 3D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 \param Szz Szz wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz){
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    
    Sxx +=sumHorizonalDerivative;
    Szz +=sumHorizonalDerivative;
    
    Syy.scale(setSurfaceZero);
    
}

