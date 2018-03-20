#pragma once

#include "../../Common/HostPrint.hpp"
#include "Derivatives.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief Derivatives namespace
        namespace Derivatives
        {

            //! \brief Class for the calculation of the Derivatives matrices for 3-D FD Simulations
            /*!
             * Calculation of derivative matrices for an equidistand grid
             *
             */
            template <typename ValueType>
            class FDTD2D : public Derivatives<ValueType>
            {
              public:
                //! Default constructor
                FDTD2D(){};

                //! Default destructor
                ~FDTD2D(){};

                FDTD2D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm);
                FDTD2D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) override;

                /* non-requiered matrixes */
                scai::lama::Matrix<ValueType> const &getDzf() const override;
                scai::lama::Matrix<ValueType> const &getDzb() const override;

              private:
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm) override;

                /* D*f: f=forward */
                using Derivatives<ValueType>::Dxf;
                using Derivatives<ValueType>::Dyf;
                /* D*b: b=backward */
                using Derivatives<ValueType>::Dxb;
                using Derivatives<ValueType>::Dyb;

                /* non-requiered matrixes */
                using Derivatives<ValueType>::Dzf;
                using Derivatives<ValueType>::Dzb;

                using Derivatives<ValueType>::DyfPressure;
                using Derivatives<ValueType>::DyfVelocity;
                using Derivatives<ValueType>::DybPressure;
                using Derivatives<ValueType>::DybVelocity;

                using Derivatives<ValueType>::useFreeSurface;
                using Derivatives<ValueType>::spatialFDorder;
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
