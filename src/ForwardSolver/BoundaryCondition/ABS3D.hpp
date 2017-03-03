#pragma once

#include "../../Acquisition/Coordinates.hpp"
#include "../../Common/HostPrint.hpp"

#include "ABS.hpp"

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

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface) override;

                void apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3, scai::lama::Vector &v4);
                void apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3, scai::lama::Vector &v4, scai::lama::Vector &v5, scai::lama::Vector &v6, scai::lama::Vector &v7, scai::lama::Vector &v8, scai::lama::Vector &v9);

              private:
                scai::lama::DenseVector<ValueType> damping; //!< Absorbing Coefficient vector
                using ABS<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
