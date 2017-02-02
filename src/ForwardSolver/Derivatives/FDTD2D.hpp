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

                FDTD2D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm);
                FDTD2D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration const &config, dmemo::CommunicatorPtr comm);

                void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration const &config, dmemo::CommunicatorPtr comm) override;

                /* non-requiered matrixes */
                lama::Matrix const &getDzf() const override;
                lama::Matrix const &getDzb() const override;

              private:
                void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm) override;

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

//! \brief Constructor to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::FDTD2D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration const &config, dmemo::CommunicatorPtr comm)
{
    init(dist, ctx, config, comm);
}

//! \brief Initialisation to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration const &config, dmemo::CommunicatorPtr comm)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");
    Derivatives<ValueType>::initializeMatrices(dist, ctx, config, comm);
}

//! \brief Constructor of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval#
 \param spatialFDorderInput FD-order of spatial derivative stencils
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::FDTD2D(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm)
{
    initializeMatrices(dist, ctx, NX, NY, NZ, DH, DT, spatialFDorderInput, comm);
}

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DH Grid spacing (equidistant)
 \param DT Temporal sampling interval
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorderInput, dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices")

    HOST_PRINT(comm, "Initialization of the matrices Dxf, Dyf, Dxb and Dyb…\n");

    // Set FD-order to class member
    spatialFDorder = spatialFDorderInput;

    /* Set FD-Coefficients */
    this->setFDCoef(spatialFDorder);

    this->calcDxf(NX, NY, NZ, dist);
    this->calcDyf(NX, NY, NZ, dist);

    HOST_PRINT(comm, "Matrix Dxf and Dyf finished.\n");

    Dxf.setContextPtr(ctx);
    Dxb.setContextPtr(ctx);
    Dyf.setContextPtr(ctx);
    Dyb.setContextPtr(ctx);

    Dxb.assignTranspose(Dxf);
    Dxb.scale(-1.0);
    Dyb.assignTranspose(Dyf);
    Dyb.scale(-1.0);

    HOST_PRINT(comm, "Matrix Dxb and Dyb finished.\n");

    Dxf.scale(lama::Scalar(DT / DH));
    Dxb.scale(lama::Scalar(DT / DH));
    Dyf.scale(lama::Scalar(DT / DH));
    Dyb.scale(lama::Scalar(DT / DH));

    HOST_PRINT(comm, "Finished with initialization of the matrices!\n");
}

//! \brief Getter method for derivative matrix Dzb
template <typename ValueType>
lama::Matrix const &KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getDzb() const
{
    COMMON_THROWEXCEPTION("There is no Dzb derivative matrix in the 2D elastic case.")
    return (Dzb);
}

//! \brief Getter method for derivative matrix Dzf
template <typename ValueType>
lama::Matrix const &KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getDzf() const
{
    COMMON_THROWEXCEPTION("There is no Dzf derivative matrix in the 2D elastic case.")
    return (Dzf);
}
