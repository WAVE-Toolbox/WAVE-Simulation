#pragma once

#include "../../Acquisition/Coordinates.hpp"
#include "../../Common/HostPrint.hpp"
#include "ABS.hpp"
#include "../../Acquisition/Coordinates.hpp"

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
            class ABS2D : public ABS<ValueType>
            {
              public:
                //! Default constructor
                ABS2D(){};

                //! Default destructor
                ~ABS2D(){};

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface) override;

                void apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3);
                void apply(scai::lama::Vector &v1, scai::lama::Vector &v2, scai::lama::Vector &v3, scai::lama::Vector &v4, scai::lama::Vector &v5);

              private:
                scai::lama::DenseVector<ValueType> damping; //!< Absorbing Coefficient vector
                using ABS<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
