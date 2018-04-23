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
            class FDTD3D : public Derivatives<ValueType>
            {
              public:
                //! Default constructor
                FDTD3D(){};

                //! Default destructor
                ~FDTD3D(){};

                FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm);
                FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) override;

              private:
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm) override;

                /* D*f: f=forward */
                using Derivatives<ValueType>::Dxf;
                using Derivatives<ValueType>::Dyf;
                using Derivatives<ValueType>::Dzf;
                /* D*b: b=backward */
                using Derivatives<ValueType>::Dxb;
                using Derivatives<ValueType>::Dyb;
                using Derivatives<ValueType>::Dzb;

                using Derivatives<ValueType>::DyfFreeSurface;
                using Derivatives<ValueType>::DybFreeSurface;

                using Derivatives<ValueType>::useFreeSurface;
                using Derivatives<ValueType>::spatialFDorder;
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
