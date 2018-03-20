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

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface) override;

                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4);
                void apply(scai::lama::Vector<ValueType> &v1, scai::lama::Vector<ValueType> &v2, scai::lama::Vector<ValueType> &v3, scai::lama::Vector<ValueType> &v4, scai::lama::Vector<ValueType> &v5, scai::lama::Vector<ValueType> &v6, scai::lama::Vector<ValueType> &v7, scai::lama::Vector<ValueType> &v8, scai::lama::Vector<ValueType> &v9);

              private:
                scai::lama::DenseVector<ValueType> damping; //!< Absorbing Coefficient vector
                using ABS<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
