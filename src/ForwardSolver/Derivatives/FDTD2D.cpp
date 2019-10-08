#include "FDTD2D.hpp"
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
KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::FDTD2D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    if (isSetup) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from init in FDTD2D.cpp : Existing Derivative Setup will be overwritten\n");
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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_ASSERT(isSetup, "call setup function before init");
    if (useSparse)
        initializeMatrices(dist, ctx, modelCoordinates, comm);
    else
        initializeMatrices(dist, ctx, modelCoordinates.getDH(), comm);

    if (useFreeSurface == 1)
        initializeFreeSurfaceMatrices(dist, ctx, modelCoordinates, comm);
}

//! \brief redistribution of all matrices
/*!
 * Only implemented for non stencil matrices. 
 * stencil matrices are only working with grid distribution. No redistribute is necessary
 * 
 \param dist Distribution of the wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::redistributeMatrices(scai::dmemo::DistributionPtr dist)
{
    SCAI_ASSERT_ERROR(Dxf.getNumRows() == 0, "redistribute isn't implemented for stencil matrices")

    DxfSparse.redistribute(dist, dist);
    DxbSparse.redistribute(dist, dist);
    DybSparse.redistribute(dist, dist);

    // if acoustic && useFreesurface==1 or if elastic && varGrid && useFreesurface==1, Dyf is not used
    if (DyfSparse.getNumRows() != 0) {
        DyfSparse.redistribute(dist, dist);
    }

    if (useVarGrid) {
        InterpolationFull.redistribute(dist, dist);
        InterpolationStaggeredX.redistribute(dist, dist);
        if (isElastic) {
            DyfStaggeredXSparse.redistribute(dist, dist);
            if (useFreeSurface != 1) { // if elastic && varGrid && useFreesurface==1, DybStaggeredXSparse is not used
                DybStaggeredXSparse.redistribute(dist, dist);
            }
        }
    }

    if (useFreeSurface == 1) {
        this->getDyfFreeSurface().redistribute(dist, dist);

        if (isElastic) {
            if (!useVarGrid) {
                this->getDybFreeSurface().redistribute(dist, dist);

            } else {
                DybStaggeredXFreeSurface.redistribute(dist, dist);
            }
        }
    }
}

template <typename ValueType>
scai::IndexType KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getNumDMatrices()
{
    IndexType numDMatrices = 4;
    if (((useVarGrid) || (useFreeSurface == 1)) && (isElastic)) {
        numDMatrices += 2;
    }
    // in acoustic modelling Dyf will be exchanged by DyfFreesurface no extra matrix is needed
    return numDMatrices;
}

template <typename ValueType>
scai::IndexType KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getNumInterpMatrices()
{
    IndexType numInterpMatrices = 0;
    if (useVarGrid) {
        // elastic: p and vx will be interpolated at the interface
        numInterpMatrices += 2;
    }

    return numInterpMatrices;
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    if (isSetup) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from estimateMemory in FDTD2D.cpp : Existing Derivative Setup will be overwritten\n");
    }
    this->setup(config);

    return (estimateMemory(dist, modelCoordinates));
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::estimateMemory(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_ASSERT(isSetup, "call setup function before estimateMemory");

    return (this->getMemoryUsage(dist, modelCoordinates, getNumDMatrices(), getNumInterpMatrices()));
}
//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("Derivatives.FDTD2D.initializeMatricesConst")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dxb and Dyb \n");

    this->calcDxf(dist);
    this->calcDyf(dist);

    HOST_PRINT(comm, "", "Matrix Dxf and Dyf finished.\n");
    Dxf.setContextPtr(ctx);
    Dxb.setContextPtr(ctx);
    Dyf.setContextPtr(ctx);
    Dyb.setContextPtr(ctx);

    Dxb = transpose(Dxf);
    Dxb.scale(-1.0);
    Dyb = transpose(Dyf);
    Dyb.scale(-1.0);

    HOST_PRINT(comm, "", "Matrix Dxb and Dyb finished.\n");
    Dxf *= this->DT / DH;
    Dxb *= this->DT / DH;
    Dyf *= this->DT / DH;
    Dyb *= this->DT / DH;

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_REGION("Derivatives.FDTD2D.initializeMatricesVar")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dxb and Dyb \n");

    this->calcDxf(modelCoordinates, dist);
    this->calcDyf(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Matrix Dxf and Dyf finished.\n");
    this->calcDxb(modelCoordinates, dist);
    this->calcDyb(modelCoordinates, dist);
    HOST_PRINT(comm, "", "Matrix Dxb and Dyb finished.\n");

    DxfSparse.setContextPtr(ctx);
    DxbSparse.setContextPtr(ctx);
    DyfSparse.setContextPtr(ctx);
    DybSparse.setContextPtr(ctx);

    DxfSparse *= this->DT;
    DxbSparse *= this->DT;
    DyfSparse *= this->DT;
    DybSparse *= this->DT;

    if ((isElastic) && (useVarGrid)) {
        this->calcDyfStaggeredX(modelCoordinates, dist);
        this->calcDybStaggeredX(modelCoordinates, dist);
        DyfStaggeredXSparse.setContextPtr(ctx);
        DybStaggeredXSparse.setContextPtr(ctx);
        DyfStaggeredXSparse *= this->DT;
        DybStaggeredXSparse *= this->DT;
    }

    if (useVarGrid) {
        this->calcInterpolationFull(modelCoordinates, dist);
        this->calcInterpolationStaggeredX(modelCoordinates, dist);
    }

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeFreeSurfaceMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    SCAI_REGION("Derivatives.FDTD2D.initializeFreeSurfaceMatrices")
    HOST_PRINT(comm, "", "Initialization of the free surface matrices  \n");
    this->calcDyfFreeSurface(modelCoordinates, dist);
    this->getDyfFreeSurface().setContextPtr(ctx);
    this->getDyfFreeSurface() *= this->DT;

    if (isElastic) {
        if (!useVarGrid) {
            this->calcDybFreeSurface(modelCoordinates, dist);
            this->getDybFreeSurface() *= this->DT;
            this->getDybFreeSurface().setContextPtr(ctx);

        } else {

            this->calcDybStaggeredXFreeSurface(modelCoordinates, dist);
            DybStaggeredXFreeSurface.setContextPtr(ctx);
            DybStaggeredXFreeSurface *= this->DT;
            DybStaggeredXSparse.purge(); // DybStaggeredX won't be used with  varGrid+FreeSurface
            DyfSparse.purge();           // DyfSparse won't be used with varGrid+FreeSurface
        }
    } else {
        //In acoustic modelling Dyf wont be used when free surface is used (see forwardsolver 2D/3D acoustic)
        this->getDyf().purge();
    }
    HOST_PRINT(comm, "", "Finished with initialization of the free surface matrices!\n");
}

//! \brief Getter method for derivative matrix Dzb
template <typename ValueType>
lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getDzb() const
{
    COMMON_THROWEXCEPTION("There is no Dzb derivative matrix in the 2D elastic case.")
    return (Dzb);
}

//! \brief Getter method for derivative matrix Dzf
template <typename ValueType>
lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getDzf() const
{
    COMMON_THROWEXCEPTION("There is no Dzf derivative matrix in the 2D elastic case.")
    return (Dzf);
}

//! \brief Getter method for combined derivative matrix
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getGraph(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    SCAI_REGION("Derivatives.FDTD2D.getGraph")
    SCAI_ASSERT(isSetup, "call setup function before init");

    if (DxbSparse.getNumRows() == 0) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Dxb,\n");
        this->calcDxb(modelCoordinates, dist);
    }

    if (DybSparse.getNumRows() == 0) {
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Dyb,\n");
        this->calcDyb(modelCoordinates, dist);
    }

    decltype(DxbSparse) temp( scai::hmemo::Context::getHostPtr() );
    temp = DxbSparse;
    temp += DybSparse;
    temp = transpose(temp);
    temp -= DxbSparse;
    DxbSparse.purge();
    temp -= DybSparse;
    DybSparse.purge();

    return (temp);
}
template class KITGPI::ForwardSolver::Derivatives::FDTD2D<float>;
template class KITGPI::ForwardSolver::Derivatives::FDTD2D<double>;
