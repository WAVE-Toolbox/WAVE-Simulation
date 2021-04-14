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

            //! \brief 3-D emem free surface
            template <typename ValueType>
            class FreeSurface3Demem : public FreeSurfaceEMEM<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface3Demem(){};

                //! Default destructor
                ~FreeSurface3Demem(){};

                void exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Szz);

              private:
                //   using FreeSurfaceEMEM<ValueType>::setSurfaceZero;
                using FreeSurfaceEMEM<ValueType>::scaleHorizontalUpdate;
                using FreeSurfaceEMEM<ValueType>::scaleVerticalUpdate;
                using FreeSurfaceEMEM<ValueType>::temp;
                using FreeSurfaceEMEM<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
