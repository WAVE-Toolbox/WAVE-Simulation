#pragma once

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
            class FreeSurfaceSH : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceSH(){};

                virtual ~FreeSurfaceSH(){};

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT) override;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
