#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../ForwardSolver/Derivatives/Derivatives.hpp"
#include "FreeSurfaceEMEM.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-Dememfree surface
            template <typename ValueType>
            class FreeSurface2Demem : public FreeSurfaceEMEM<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Demem(){};

                //! Default destructor
                ~FreeSurface2Demem(){};
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
