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
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{

    if (useSparse)
        initializeMatrices(dist, ctx, modelCoordinates, comm);
    else
        initializeMatrices(dist, ctx, modelCoordinates.getDH(), comm);

    if (useFreeSurface == 1)
        initializeFreeSurfaceMatrices(dist, ctx, modelCoordinates, comm);
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

    if (useVarGrid) {
        InterpolationFull.redistribute(dist, dist);
        if (isElastic) {
            DyfStaggeredXSparse.redistribute(dist, dist);
            DybStaggeredXSparse.redistribute(dist, dist);
            InterpolationStaggeredX.redistribute(dist, dist);
        }
    }

    if (useFreeSurface == 1) {
        DyfFreeSurface.redistribute(dist, dist);

        if (isElastic) {
            if (!useVarGrid) {
                DybFreeSurface.redistribute(dist, dist);

            } else {
                DybStaggeredXFreeSurface.redistribute(dist, dist);
            }
        }
    }
}

template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::FDTD2D<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    IndexType NumDMatrices = 4;
    IndexType NumInterpMatrices = 0;
    if (useVarGrid) {
        NumInterpMatrices += 1;
        if (isElastic) {
            NumDMatrices += 2;
            NumInterpMatrices += 1;
        }
    }

    if (useFreeSurface == 1) {
        NumDMatrices += 2;
    }

    if (useSparse) {
        HOST_PRINT(dist->getCommunicatorPtr(), "Total size of derivative matrices per core " << this->getMemorySparseMatrix(dist, modelCoordinates) / 1024 / 1024 * NumDMatrices / dist->getNumPartitions() << " MB"
                                                                                             << "\n",
                   "Memory of derivative Matrix: " << this->getMemorySparseMatrix(dist, modelCoordinates) / 1024 / 1024 << " MB, "
                                                   << " partitions : " << dist->getNumPartitions() << ", Num matrices : " << NumDMatrices << "\n");
    } else {
        HOST_PRINT(dist->getCommunicatorPtr(), "Total size of derivative matrices per core " << this->getMemorySparseMatrix(dist) / 1024 / 1024 * NumDMatrices / dist->getNumPartitions() << " MB"
                                                                                             << "\n",
                   "Memory of derivative Matrix: " << this->getMemoryStencilMatrix(dist) / 1024 / 1024 << " MB, "
                                                   << " partitions : " << dist->getNumPartitions() << ", Num matrices : " << NumDMatrices << "\n");
    }

    if (useVarGrid) {
        HOST_PRINT(dist->getCommunicatorPtr(), "Total size of interpolation matrices per core " << this->getMemoryInterpolationMatrix(dist) / 1024 / 1024 * NumInterpMatrices / dist->getNumPartitions() << " MB"
                                                                                                << "\n",
                   "Memory of interpolation Matrix: " << this->getMemoryInterpolationMatrix(dist) / 1024 / 1024 << " MB, "
                                                      << " partitions : " << dist->getNumPartitions() << ", Num matrices : " << NumInterpMatrices << "\n");
    }
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

    SCAI_REGION("initializeMatrices")

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
        if (isElastic) {
            this->calcInterpolationStaggeredX(modelCoordinates, dist);
        }
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
    this->calcDyfFreeSurface(modelCoordinates, dist);
    DyfFreeSurface.setContextPtr(ctx);
    DyfFreeSurface *= this->DT;

    if (isElastic) {
        if (!useVarGrid) {
            this->calcDybFreeSurface(modelCoordinates, dist);
            DybFreeSurface *= this->DT;
            DybFreeSurface.setContextPtr(ctx);

        } else {

            this->calcDybStaggeredXFreeSurface(modelCoordinates, dist);
            DybStaggeredXFreeSurface.setContextPtr(ctx);
            DybStaggeredXFreeSurface *= this->DT;
        }
    }
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
