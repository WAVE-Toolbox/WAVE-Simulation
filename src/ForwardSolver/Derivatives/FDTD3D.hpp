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

                FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, scai::IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm);
                FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm);

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm, std::vector<scai::IndexType> &FDorder) override;

                virtual void redistributeMatrices(scai::dmemo::DistributionPtr dist) override;

                scai::lama::CSRSparseMatrix<ValueType> getCombinedMatrix() override;

              private:
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm) override;
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, scai::dmemo::CommunicatorPtr comm) override;

                /* D*f: f=forward */
                using Derivatives<ValueType>::Dxf;
                using Derivatives<ValueType>::Dyf;
                using Derivatives<ValueType>::Dzf;
                using Derivatives<ValueType>::DxfSparse;
                using Derivatives<ValueType>::DyfSparse;
                using Derivatives<ValueType>::DzfSparse;
                /* D*b: b=backward */
                using Derivatives<ValueType>::Dxb;
                using Derivatives<ValueType>::Dyb;
                using Derivatives<ValueType>::Dzb;
                using Derivatives<ValueType>::DxbSparse;
                using Derivatives<ValueType>::DybSparse;
                using Derivatives<ValueType>::DzbSparse;

                using Derivatives<ValueType>::DyfStaggeredXSparse;
                using Derivatives<ValueType>::DybStaggeredXSparse;
                using Derivatives<ValueType>::DyfStaggeredZSparse;
                using Derivatives<ValueType>::DybStaggeredZSparse;
                using Derivatives<ValueType>::DyfStaggeredXZSparse;
                using Derivatives<ValueType>::DybStaggeredXZSparse;

                using Derivatives<ValueType>::DyfFreeSurface;
                using Derivatives<ValueType>::DybFreeSurface;

                using Derivatives<ValueType>::useFreeSurface;
                using Derivatives<ValueType>::useSparse;
                using Derivatives<ValueType>::useVarFDorder;
                using Derivatives<ValueType>::useVarGrid;

                using Derivatives<ValueType>::InterpolationFull;
                using Derivatives<ValueType>::InterpolationStaggeredX;
                using Derivatives<ValueType>::InterpolationStaggeredZ;
                using Derivatives<ValueType>::InterpolationStaggeredXZ;
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
