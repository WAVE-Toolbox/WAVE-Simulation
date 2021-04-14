#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../ForwardSolver/Derivatives/Derivatives.hpp"
#include "FreeSurfaceTMEM.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D tmem free surface
            template <typename ValueType>
            class FreeSurface2Dtmem : public FreeSurfaceTMEM<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Dtmem(){};

                //! Default destructor
                ~FreeSurface2Dtmem(){};

                void exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx);

              private:
                using FreeSurfaceTMEM<ValueType>::scaleHorizontalUpdate;
                using FreeSurfaceTMEM<ValueType>::scaleVerticalUpdate;
                using FreeSurfaceTMEM<ValueType>::temp;
                using FreeSurfaceTMEM<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
