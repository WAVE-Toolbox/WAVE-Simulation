#pragma once

#include "../../Common/HostPrint.hpp"
#include "ABS.hpp"

using namespace scai;

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Class for the calculation of the Absorbing Coefficient matrix for 3-D FD Simulations
            /*!
             * Calculation of the absorbing coefficient matrix for an equidistand grid
             *
             */
            template <typename ValueType>
            class ABSEM2D : public ABSEM<ValueType>
            {
              public:
                //! Default constructor
                ABSEM2D(){};

                //! Default destructor
                ~ABSEM2D(){};

                ValueType estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface) override;

                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3);
              private:
                typedef typename ABSEM<ValueType>::VectorType VectorType;
                VectorType damping; //!< Absorbing Coefficient DenseVector. damping=1.0 in the interior and  damping < 1.0 inside the boundary frame.
                using ABSEM<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
