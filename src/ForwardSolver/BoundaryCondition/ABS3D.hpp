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
            class ABS3D : public ABS<ValueType>
            {
              public:
                //! Default constructor
                ABS3D(){};

                //! Default destructor
                ~ABS3D(){};

                ValueType estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType BoundaryWidth, ValueType DampingCoeff, scai::IndexType useFreeSurface) override;

                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4);
                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4, scai::lama::Vector<ValueType> &v5, scai::lama::Vector<ValueType> &v6);
                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4, scai::lama::Vector<ValueType> &v5, scai::lama::Vector<ValueType> &v6, scai::lama::Vector<ValueType> &v7, scai::lama::Vector<ValueType> &v8, scai::lama::Vector<ValueType> &v9);

              private:
                typedef typename ABS<ValueType>::VectorType VectorType;
                VectorType damping; //!< Absorbing Coefficient DenseVector. damping=1.0 in the interior and  damping < 1.0 inside the boundary frame.
                using ABS<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
