#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurfaceSH.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D sh free surface
            template <typename ValueType>
            class FreeSurface2Dsh : public FreeSurfaceSH<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Dsh(){};

                //! Default destructor
                ~FreeSurface2Dsh(){};
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
