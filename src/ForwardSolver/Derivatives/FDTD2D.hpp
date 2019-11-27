#pragma once

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

                FDTD2D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm);

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;

                void redistributeMatrices(scai::dmemo::DistributionPtr dist) override;

                ValueType estimateMemory(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;
                ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                /* non-requiered matrixes */
                scai::lama::Matrix<ValueType> const &getDzf() const override;
                scai::lama::Matrix<ValueType> const &getDzb() const override;

                scai::lama::CSRSparseMatrix<ValueType> getGraph(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

              private:
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, scai::dmemo::CommunicatorPtr comm) override;
                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
                void initializeFreeSurfaceMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;

                scai::IndexType getNumDMatrices();
                scai::IndexType getNumInterpMatrices();

                /* D*f: f=forward */
                using Derivatives<ValueType>::Dxf;
                using Derivatives<ValueType>::Dyf;
                using Derivatives<ValueType>::DxfSparse;
                using Derivatives<ValueType>::DyfSparse;

                /* D*b: b=backward */
                using Derivatives<ValueType>::Dxb;
                using Derivatives<ValueType>::Dyb;
                using Derivatives<ValueType>::DxbSparse;
                using Derivatives<ValueType>::DybSparse;

                /* non-requiered matrixes */
                using Derivatives<ValueType>::Dzf;
                using Derivatives<ValueType>::Dzb;
                using Derivatives<ValueType>::DzfSparse;
                using Derivatives<ValueType>::DzbSparse;

                using Derivatives<ValueType>::DyfStaggeredXSparse;
                using Derivatives<ValueType>::DybStaggeredXSparse;
                /*
                using Derivatives<ValueType>::DyfFreeSurfaceSparse;
                using Derivatives<ValueType>::DybFreeSurfaceSparse;*/
                using Derivatives<ValueType>::DybStaggeredXFreeSurface;

                using Derivatives<ValueType>::useFreeSurface;
                using Derivatives<ValueType>::useStencilMatrix;
                using Derivatives<ValueType>::useVarGrid;
                using Derivatives<ValueType>::isElastic;

                using Derivatives<ValueType>::useVarFDorder;
                using Derivatives<ValueType>::isSetup;

                using Derivatives<ValueType>::InterpolationFull;
                using Derivatives<ValueType>::InterpolationStaggeredX;
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
