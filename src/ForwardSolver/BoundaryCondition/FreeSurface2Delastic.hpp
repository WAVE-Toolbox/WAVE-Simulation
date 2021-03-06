#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurfaceElastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D elastic free surface
            template <typename ValueType>
            class FreeSurface2Delastic : public FreeSurfaceElastic<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Delastic(){};

                //! Default destructor
                ~FreeSurface2Delastic(){};

                void exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx);

              private:
                using FreeSurfaceElastic<ValueType>::scaleHorizontalUpdate;
                using FreeSurfaceElastic<ValueType>::scaleVerticalUpdate;
                using FreeSurfaceElastic<ValueType>::temp;
                using FreeSurfaceElastic<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
