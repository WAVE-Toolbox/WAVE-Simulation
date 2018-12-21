#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>
#include "../../Acquisition/Coordinates.hpp"

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
                 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
                 \param BoundaryWidth Width of damping boundary
                 \param DampingCoeff Damping coefficient
                 \param useFreeSurface Indicator which free surface is in use
                 */
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates const &modelCoordinates, scai::IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface) = 0;

              protected:
                // For the ABS Boundaries Sparse Vectors and Dense Vectors can be declared. The code will run without any further changes.
                typedef typename scai::lama::DenseVector<ValueType> VectorType; //!< Define Vector Type as Dense vector. For big models switch to SparseVector
            };
        } /* end namespace BoundaryCondition  */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
