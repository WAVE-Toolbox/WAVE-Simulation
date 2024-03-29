#include "FDTD3D.hpp"
#include <scai/common/Settings.hpp>

using namespace scai;

//! \brief Constructor to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    this->setup(config);
    init(dist, ctx, modelCoordinates, comm);
}

//! \brief Initialisation to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (isSetup) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from init in FDTD3D.cpp : Existing Derivative Setup will be overwritten\n");
    }
    this->setup(config);

    init(dist, ctx, modelCoordinates, comm);
}

//! \brief Initialisation to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_ASSERT(isSetup, "call setup function before init");
    if (useStencilMatrix) {
        initializeMatrices(dist, ctx, modelCoordinates.getDH(), comm);
    } else {
        initializeMatrices(dist, ctx, modelCoordinates, comm);
    }
    if (useFreeSurface == 1)
        initializeFreeSurfaceMatrices(dist, ctx, modelCoordinates, comm);
}

//! \brief redistribution of all matrices
/*!
 *
 \param dist Distribution of the wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::redistributeMatrices(scai::dmemo::DistributionPtr dist)
{
    HOST_PRINT(dist->getCommunicatorPtr(), "", "started redistribution of the derivative matrices.\n");

    SCAI_ASSERT_ERROR(Dxf.getNumRows() == 0, "redistribute isn't implemented for stencil matrices")

    DxfSparse.redistribute(dist, dist);

    // if acoustic && useFreesurface==1 or if elastic && varGrid && useFreesurface==1, Dyf is not used
    if (DyfSparse.getNumRows() != 0)
        DyfSparse.redistribute(dist, dist);

    DzfSparse.redistribute(dist, dist);
    DxbSparse.redistribute(dist, dist);
    DybSparse.redistribute(dist, dist);
    DzbSparse.redistribute(dist, dist);

    if (useVarGrid) {
        InterpolationStaggeredX.redistribute(dist, dist);
        InterpolationStaggeredZ.redistribute(dist, dist);
        InterpolationFull.redistribute(dist, dist);
        if (isElastic) {
            DyfStaggeredXSparse.redistribute(dist, dist);
            DyfStaggeredZSparse.redistribute(dist, dist);

            if (useFreeSurface != 1) { // if elastic && varGrid && useFreesurface==1, DybStaggeredXSparse nad DybStaggeredZSparse are not used
                DybStaggeredXSparse.redistribute(dist, dist);
                DybStaggeredZSparse.redistribute(dist, dist);
            }

            InterpolationStaggeredXZ.redistribute(dist, dist);
        }
    }

    if (useFreeSurface == 1) {
        this->getDyfFreeSurface().redistribute(dist, dist);

        if (isElastic) {
            if (!useVarGrid) {
                this->getDybFreeSurface().redistribute(dist, dist);

            } else {
                DybStaggeredXFreeSurface.redistribute(dist, dist);
                DybStaggeredZFreeSurface.redistribute(dist, dist);
            }
        }
    }
    HOST_PRINT(dist->getCommunicatorPtr(), "", "finished redistribution of the derivative matrices.\n");
}

template <typename ValueType>
scai::IndexType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::getNumDMatrices()
{
    IndexType numDMatrices = 6;

    if ((useVarGrid) && (isElastic)) {
        /* + DyfStaggeredX + DybStaggeredX + DyfStaggeredZ + DybStaggeredZ
            or with free surface:
            +DyfFreeSurface - Dyf +  DyfStaggeredX + DyfStaggeredZ + DybStaggeredXFreeSurface + DybStaggeredZFreeSurface
            */
        numDMatrices += 4;
    }

    if ((useFreeSurface == 1) && (isElastic) && (!useVarGrid)) {
        // + DyfFreeSurface and DybFreeSurface/StaggeredX -
        numDMatrices += 2;
    }
    return numDMatrices;
}

template <typename ValueType>
scai::IndexType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::getNumInterpMatrices()
{
    IndexType numInterpMatrices = 0;
    if (useVarGrid) {
        // + interpolation for (a^  Sxx,Syy,Szz,P),Vx,Vz,
        numInterpMatrices += 3;
        if (isElastic) {
            // + Interpolation for Sxz
            numInterpMatrices += 1;
        }
    }

    return numInterpMatrices;
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    if (isSetup) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from estimateMemory FDTD3D.cpp : Existing Derivative Setup will be overwritten\n");
    }
    this->setup(config);

    return (estimateMemory(dist, modelCoordinates));
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::estimateMemory(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{

    SCAI_ASSERT(isSetup, "call setup function before estimateMemory");
    return (this->getMemoryUsage(dist, modelCoordinates, getNumDMatrices(), getNumInterpMatrices()));
}

//! \brief Initialization of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_REGION("Derivatives.FDTD3D.initializeMatricesConst")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dzf, Dxb, Dyb, Dzb \n");

    this->calcDxf(dist);
    this->calcDzf(dist);
    this->calcDyf(dist);

    HOST_PRINT(comm, "", "Matrix Dxf, Dyf and Dzf finished.\n");

    Dxf.setContextPtr(ctx);
    Dzf.setContextPtr(ctx);
    Dxb.setContextPtr(ctx);
    Dzb.setContextPtr(ctx);
    Dyf.setContextPtr(ctx);
    Dyb.setContextPtr(ctx);

    Dxb.assignTranspose(Dxf);
    Dxb *= -1.0;
    Dzb.assignTranspose(Dzf);
    Dzb *= -1.0;
    Dyb.assignTranspose(Dyf);
    Dyb *= -1.0;

    HOST_PRINT(comm, "", "Matrix Dxb, Dyb and Dzb finished.\n");

    Dxf.scale(this->DT / DH);
    Dzf.scale(this->DT / DH);
    Dxb.scale(this->DT / DH);
    Dzb.scale(this->DT / DH);
    Dyf.scale(this->DT / DH);
    Dyb.scale(this->DT / DH);

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Initialization of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_REGION("Derivatives.FDTD3D.initializeMatricesVar")

    HOST_PRINT(comm, "Initialization of the matrices: \n");
    HOST_PRINT(comm, "", "Dxf," << std::flush);
    this->calcDxf(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Dzf," << std::flush);
    this->calcDzf(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Dxb," << std::flush);
    this->calcDxb(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Dyb," << std::flush);
    this->calcDyb(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Dzb," << std::flush);
    this->calcDzb(modelCoordinates, dist);

    DxfSparse.setContextPtr(ctx);
    DxbSparse.setContextPtr(ctx);
    DybSparse.setContextPtr(ctx);
    DzfSparse.setContextPtr(ctx);
    DzbSparse.setContextPtr(ctx);

    DxfSparse.scale(this->DT);
    DxbSparse.scale(this->DT);
    DybSparse.scale(this->DT);
    DzfSparse.scale(this->DT);
    DzbSparse.scale(this->DT);

    if (!((useFreeSurface == 1) && (useVarGrid))) {
        HOST_PRINT(comm, "", "Dyf," << std::flush);
        this->calcDyf(modelCoordinates, dist);
        DyfSparse.setContextPtr(ctx);
        DyfSparse.scale(this->DT);
    }

    if ((isElastic) && (useVarGrid)) {
        HOST_PRINT(comm, "", "DyfStaggeredX," << std::flush);
        this->calcDyfStaggeredX(modelCoordinates, dist);

        HOST_PRINT(comm, "", "DyfStaggeredZ," << std::flush);
        this->calcDyfStaggeredZ(modelCoordinates, dist);

        DyfStaggeredXSparse.setContextPtr(ctx);
        DyfStaggeredZSparse.setContextPtr(ctx);

        DyfStaggeredXSparse *= this->DT;
        DyfStaggeredZSparse *= this->DT;

        if (useFreeSurface != 1) {
            HOST_PRINT(comm, "", "DybStaggeredX," << std::flush);
            this->calcDybStaggeredX(modelCoordinates, dist);
            HOST_PRINT(comm, "", "DybStaggeredZ," << std::flush);
            this->calcDybStaggeredZ(modelCoordinates, dist);

            DybStaggeredXSparse.setContextPtr(ctx);
            DybStaggeredZSparse.setContextPtr(ctx);
            DybStaggeredXSparse *= this->DT;
            DybStaggeredZSparse *= this->DT;
        }
    }

    if (useVarGrid) {
        this->calcInterpolationFull(modelCoordinates, dist);
        this->calcInterpolationStaggeredX(modelCoordinates, dist);
        this->calcInterpolationStaggeredZ(modelCoordinates, dist);
        if (isElastic) {
            this->calcInterpolationStaggeredXZ(modelCoordinates, dist);
        }
    }

    //HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Initialization of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeFreeSurfaceMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_REGION("Derivatives.FDTD3D.initializeFreeSurfaceMatrices")

    HOST_PRINT(comm, "", "DyfFreeSurface," << std::flush);
    this->calcDyfFreeSurface(modelCoordinates, dist);
    this->getDyfFreeSurface().setContextPtr(ctx);
    this->getDyfFreeSurface() *= this->DT;

    if (isElastic) {
        if (!useVarGrid) {
            HOST_PRINT(comm, "", "DybFreeSurface," << std::flush);
            this->calcDybFreeSurface(modelCoordinates, dist);
            this->getDybFreeSurface().setContextPtr(ctx);
            this->getDybFreeSurface() *= this->DT;

        } else {
            //DyfSparse.purge(); // DyfSparse won't be used with varGrid+FreeSurface
            HOST_PRINT(comm, "", "DybStaggeredXFreeSurface," << std::flush);
            this->calcDybStaggeredXFreeSurface(modelCoordinates, dist);
            HOST_PRINT(comm, "", "DybStaggeredZFreeSurface," << std::flush);
            this->calcDybStaggeredZFreeSurface(modelCoordinates, dist);
            DybStaggeredXFreeSurface.setContextPtr(ctx);
            DybStaggeredZFreeSurface.setContextPtr(ctx);
            DybStaggeredXFreeSurface *= this->DT;
            DybStaggeredZFreeSurface *= this->DT;
        }
    } else {
        //In acoustic modelling Dyf wont be used when free surface is used (see forwardsolver 2D/3D acoustic)
        this->getDyf().purge();
    }
}

//! \brief Getter method for combined derivative matrix
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::getGraph(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_REGION("Derivatives.FDTD3D.getGraph")
    SCAI_ASSERT(isSetup, "call setup function before init");

    if (DxbSparse.getNumRows() == 0) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Dxb," << std::flush);
        this->calcDxb(modelCoordinates, dist);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Dyb," << std::flush);
        this->calcDyb(modelCoordinates, dist);
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Dzb," << std::flush);
        this->calcDzb(modelCoordinates, dist);
    }

    decltype(DxbSparse) temp(scai::hmemo::Context::getHostPtr());
    temp = DxbSparse;
    temp += DybSparse;
    temp += DzbSparse;
    temp = transpose(temp);
    temp -= DxbSparse;
    DxbSparse.purge();
    temp -= DybSparse;
    DybSparse.purge();
    temp -= DzbSparse;
    DzbSparse.purge();

    return (temp);
}

template class KITGPI::ForwardSolver::Derivatives::FDTD3D<double>;
template class KITGPI::ForwardSolver::Derivatives::FDTD3D<float>;
