#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../Modelparameter/Modelparameter.hpp"
#include "../../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../../ForwardSolver/Derivatives/FDTD3D.hpp"
#include "FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D emem free surface
            template <typename ValueType>
            class FreeSurfaceEMEM : public FreeSurfaceEM<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceEMEM(){};

                virtual ~FreeSurfaceEMEM() = 0;

                void init(scai::dmemo::DistributionPtr dist, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType /*DT*/) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model);

              protected:
                using FreeSurfaceEM<ValueType>::active;
                using FreeSurfaceEM<ValueType>::setZeroFreeSurface;
                scai::lama::SparseVector<ValueType> selectFreeSurface; //!< //!< Vector, which sets everything besides the free surface to zero
                scai::lama::SparseVector<ValueType> temp;
                scai::lama::SparseVector<ValueType> scaleHorizontalUpdate; //!< Vector, which scales the horizontal update on the free surface in order to exchange the horizontal main stresses.
                scai::lama::SparseVector<ValueType> scaleVerticalUpdate;   //!< Vector, which scales the vertical updateon the free surface in order to exchange the horizontal main stresses.
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
