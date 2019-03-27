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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");

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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm, std::vector<IndexType> &FDorder)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");
    useSparse = true;
    useVarFDorder = true;

    this->setFDCoef();

    this->setFDOrder(config.get<std::string>("spatialFDorderFilename"));

    initializeMatrices(dist, ctx, modelCoordinates, config.get<ValueType>("DT"), comm);
}

//! \brief redistribution of all matrices
/*!
 *
 \param dist Distribution of the wavefield
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::redistributeMatrices(scai::dmemo::DistributionPtr dist)
{
    DxfSparse.redistribute(dist, dist);
    DyfSparse.redistribute(dist, dist);
    DxbSparse.redistribute(dist, dist);
    DybSparse.redistribute(dist, dist);
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
KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::FDTD2D(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, IndexType spatialFDorderInput, scai::dmemo::CommunicatorPtr comm)
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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices")

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
    Dxf *= DT / DH;
    Dxb *= DT / DH;
    Dyf *= DT / DH;
    Dyb *= DT / DH;

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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dxb and Dyb \n");

    this->calcDxf(modelCoordinates, dist);
    this->calcDyf(modelCoordinates, dist);
    this->calcDxb(modelCoordinates, dist);
    this->calcDyb(modelCoordinates, dist);

    HOST_PRINT(comm, "", "Matrix Dxf and Dyf finished.\n");
    DxfSparse.setContextPtr(ctx);
    DxbSparse.setContextPtr(ctx);
    DyfSparse.setContextPtr(ctx);
    DybSparse.setContextPtr(ctx);

    HOST_PRINT(comm, "", "Matrix Dxb and Dyb finished.\n");
    DxfSparse *= DT;
    DxbSparse *= DT;
    DyfSparse *= DT;
    DybSparse *= DT;

    if (modelCoordinates.isVariable())
        this->calcInterpolationP(modelCoordinates, dist);

    HOST_PRINT(comm, "", "Finished with initialization of the matrices!\n");
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
scai::lama::CSRSparseMatrix<ValueType> KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::getCombinedMatrix()
{
    auto temp = DxfSparse;
    temp += DyfSparse;
    temp -= DxbSparse;
    temp -= DybSparse;

    return (temp);
}
template class KITGPI::ForwardSolver::Derivatives::FDTD2D<float>;
template class KITGPI::ForwardSolver::Derivatives::FDTD2D<double>;
