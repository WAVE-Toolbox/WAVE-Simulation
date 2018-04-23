#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurfaceVisco.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D visco free surface
            template <typename ValueType>
            class FreeSurface2Dvisco : public FreeSurfaceVisco<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Dvisco(){};

                //! Default destructor
                ~FreeSurface2Dvisco(){};

                void exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Rxx, ValueType DThalf);

              private:
                using FreeSurfaceVisco<ValueType>::scaleStressHorizontalUpdate;     //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                using FreeSurfaceVisco<ValueType>::scaleStressVerticalUpdate;       //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                using FreeSurfaceVisco<ValueType>::scaleRelaxationHorizontalUpdate; //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                using FreeSurfaceVisco<ValueType>::scaleRelaxationVerticalUpdate;   //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                using FreeSurfaceVisco<ValueType>::selectHorizontalUpdate;          //!< Vector, which sets everything besides the free surface to zero
                using FreeSurfaceVisco<ValueType>::temp;                            //!< temporary Sparse Vector (lives as long as the forward solver)
                using FreeSurfaceVisco<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
