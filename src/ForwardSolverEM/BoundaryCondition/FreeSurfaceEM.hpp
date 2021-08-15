#pragma once

#include "../../ForwardSolver/BoundaryCondition/FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief Abstract free surface class
            template <typename ValueType>
            class FreeSurfaceEM : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceEM(){};

                //! Default destructor
                virtual ~FreeSurfaceEM(){};

                virtual void init(scai::dmemo::DistributionPtr dist, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT) override;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
