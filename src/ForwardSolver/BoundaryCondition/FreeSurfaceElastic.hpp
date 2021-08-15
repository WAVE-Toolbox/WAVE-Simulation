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

            //! \brief 3-D elastic free surface
            template <typename ValueType>
            class FreeSurfaceElastic : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceElastic(){};

                virtual ~FreeSurfaceElastic(){};

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType /*DT*/) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model);

              protected:
                using FreeSurface<ValueType>::active;
                using FreeSurface<ValueType>::setZeroFreeSurface;
                scai::lama::SparseVector<ValueType> selectFreeSurface; //!< //!< Vector, which sets everything besides the free surface to zero
                scai::lama::SparseVector<ValueType> temp;
                scai::lama::SparseVector<ValueType> scaleHorizontalUpdate; //!< Vector, which scales the horizontal update on the free surface in order to exchange the horizontal main stresses.
                scai::lama::SparseVector<ValueType> scaleVerticalUpdate;   //!< Vector, which scales the vertical updateon the free surface in order to exchange the horizontal main stresses.
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
