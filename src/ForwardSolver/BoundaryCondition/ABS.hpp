#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

using namespace scai;

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Abstract class for the calculation of the Absorbing Boundary condition
            template <typename ValueType>
            class ABS
            {
              public:
                //! \brief Default constructor
                ABS()
                    : active(false){};

                //! \brief Default destructor
                ~ABS(){};

                bool active; //!< Bool if ABS method is active

                //! \brief Initializsation of the Absorbing-Coefficient matrix
                /*!
                 *
                 \param dist Distribution of the wavefield
                 \param ctx Context
                 \param NX Total number of grid points in X
                 \param NY Total number of grid points in Y
                 \param NZ Total number of grid points in Z
                 \param BoundaryWidth Width of damping boundary
                 \param DampingCoeff Damping coefficient
                 \param useFreeSurface Bool if free surface is in use
                 */
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface) = 0;
            };
        } /* end namespace BoundaryCondition  */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
