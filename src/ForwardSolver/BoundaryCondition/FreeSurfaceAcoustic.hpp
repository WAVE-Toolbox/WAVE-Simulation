#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../Modelparameter/Modelparameter.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "../Derivatives/FDTD3D.hpp"
#include "FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D acoustic free surface
            template <typename ValueType>
            class FreeSurfaceAcoustic : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceAcoustic(){};

                virtual ~FreeSurfaceAcoustic() = 0;

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT) override;

              protected:
                using FreeSurface<ValueType>::active;
                using FreeSurface<ValueType>::setZeroFreeSurface;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
