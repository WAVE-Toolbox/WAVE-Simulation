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
    init(dist, ctx, config, modelCoordinates, comm);
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
    useFreeSurface = config.get<IndexType>("FreeSurface");
    useVarGrid = modelCoordinates.isVariable();

    this->setFDCoef();

    auto useVariableGrid = config.get<bool>("useVariableGrid");
    auto useGraphPartitioning = config.get<bool>("useGraphPartitioning");

    useFreeSurface = config.get<IndexType>("FreeSurface");
    if ((useVariableGrid) || (useGraphPartitioning)) {
        useSparse = true;
    }

    if ((useVariableGrid || useGraphPartitioning) && config.get<bool>("useVariableFDoperators")) {
        useVarFDorder = true;
        this->setFDOrder(config.get<std::string>("spatialFDorderFilename"));
    } else {
        // Set FD-order to class member
        this->setFDOrder(config.get<IndexType>("spatialFDorder"));
    }

    if (useSparse)
        initializeMatrices(dist, ctx, modelCoordinates, config.get<ValueType>("DT"), comm);
    else
        initializeMatrices(dist, ctx, modelCoordinates.getDH(), config.get<ValueType>("DT"), comm);
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
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm, std::vector<IndexType> &FDorder)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");
    useVarGrid = modelCoordinates.isVariable();
    useSparse = true;
    useVarFDorder = true;

    this->setFDCoef();
    this->setFDOrder(FDorder);

    initializeMatrices(dist, ctx, modelCoordinates, config.get<ValueType>("DT"), comm);
}

//! \brief redistribution of all matrices
/*!
 *
 \param dist Distribution of the wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::redistributeMatrices(scai::dmemo::DistributionPtr dist)
{
    DxfSparse.redistribute(dist, dist);
    DyfSparse.redistribute(dist, dist);
    DzfSparse.redistribute(dist, dist);
    DxbSparse.redistribute(dist, dist);
    DybSparse.redistribute(dist, dist);
    DzbSparse.redistribute(dist, dist);
    if (useVarGrid)
        InterpolationP.redistribute(dist, dist);
}

//! \brief Constructor of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal sampling interval#
 \param spatialFDorderInput FD-order of spatial derivative stencils
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::FDTD3D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm)
{
    initializeMatrices(dist, ctx, modelCoordinates, DT, comm);
}

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param DT Temporal sampling interval
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices");

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

    lama::SyncKind syncKind = lama::SyncKind::SYNCHRONOUS;

    // by default do matrix-vector operations synchronously, but can be set via environment
    bool isSet;

    if (common::Settings::getEnvironment(isSet, "SCAI_ASYNCHRONOUS")) {
        if (isSet) {
            syncKind = lama::SyncKind::ASYNC_COMM;
        }
    }

    Dxf.setCommunicationKind(syncKind);
    Dzf.setCommunicationKind(syncKind);
    Dxb.setCommunicationKind(syncKind);
    Dzb.setCommunicationKind(syncKind);
    Dyf.setCommunicationKind(syncKind);
    Dyb.setCommunicationKind(syncKind);

    Dxb.assignTranspose(Dxf);
    Dxb *= -1.0;
    Dzb.assignTranspose(Dzf);
    Dzb *= -1.0;
    Dyb.assignTranspose(Dyf);
    Dyb *= -1.0;

    HOST_PRINT(comm, "", "Matrix Dxb, Dyb and Dzb finished.\n");

    Dxf.scale(DT / DH);
    Dzf.scale(DT / DH);
    Dxb.scale(DT / DH);
    Dzb.scale(DT / DH);
    Dyf.scale(DT / DH);
    Dyb.scale(DT / DH);

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Initializsation of the derivative matrices
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal sampling interval
 \param spatialFDorderInput FD-order of spatial stencils
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dzf, Dxb, Dyb, Dzb\n");

    this->calcDxf(modelCoordinates, dist);
    this->calcDyf(modelCoordinates, dist);
    this->calcDzf(modelCoordinates, dist);
    this->calcDxb(modelCoordinates, dist);
    this->calcDyb(modelCoordinates, dist);
    this->calcDzb(modelCoordinates, dist);

    HOST_PRINT(comm, "", "Matrix Dxf, Dyf and Dzf finished.\n");

    DxfSparse.setContextPtr(ctx);
    DxbSparse.setContextPtr(ctx);
    DyfSparse.setContextPtr(ctx);
    DybSparse.setContextPtr(ctx);
    DzfSparse.setContextPtr(ctx);
    DzbSparse.setContextPtr(ctx);

    lama::SyncKind syncKind = lama::SyncKind::SYNCHRONOUS;

    // by default do matrix-vector operations synchronously, but can be set via environment
    bool isSet;

    if (common::Settings::getEnvironment(isSet, "SCAI_ASYNCHRONOUS")) {
        if (isSet) {
            syncKind = lama::SyncKind::ASYNC_COMM;
        }
    }

    DxfSparse.setCommunicationKind(syncKind);
    DxbSparse.setCommunicationKind(syncKind);
    DyfSparse.setCommunicationKind(syncKind);
    DybSparse.setCommunicationKind(syncKind);
    DzfSparse.setCommunicationKind(syncKind);
    DzbSparse.setCommunicationKind(syncKind);

    HOST_PRINT(comm, "", "Matrix Dxb, Dyb and Dzb finished.\n");

    DxfSparse.scale(DT);
    DxbSparse.scale(DT);
    DyfSparse.scale(DT);
    DybSparse.scale(DT);
    DzfSparse.scale(DT);
    DzbSparse.scale(DT);

    if (modelCoordinates.isVariable())
        this->calcInterpolationP(modelCoordinates, dist);

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
}

//! \brief Getter method for combined derivative matrix
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::getCombinedMatrix()
{
    auto temp = DxfSparse;
    temp += DyfSparse;
    temp += DzfSparse;
    temp -= DxbSparse;
    temp -= DybSparse;
    temp -= DzbSparse;
    return (temp);
}

template class KITGPI::ForwardSolver::Derivatives::FDTD3D<double>;
template class KITGPI::ForwardSolver::Derivatives::FDTD3D<float>;
