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
    if(isSetup){
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from init in FDTD3D.cpp : Existing Derivative Setup will be overwritten\n");
    }
    this->setup(config);
    
    init(dist,ctx,modelCoordinates,comm);
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
    if (useSparse) {
        initializeMatrices(dist, ctx, modelCoordinates, comm);
    } else {
        initializeMatrices(dist, ctx, modelCoordinates.getDH(), comm);
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
    
    SCAI_ASSERT_ERROR(Dxf.getNumRows()==0, "redistribute isn't implemented for stencil matrices")
   
    DxfSparse.redistribute(dist, dist);
    
    // if acoustic && useFreesurface==1 or if elastic && varGrid && useFreesurface==1, Dyf is not used
    if(DyfSparse.getNumRows()!=0)
    DyfSparse.redistribute(dist, dist);
    
    DzfSparse.redistribute(dist, dist);
    DxbSparse.redistribute(dist, dist);
    DybSparse.redistribute(dist, dist);
    DzbSparse.redistribute(dist, dist);

    if (useVarGrid) {
        if (isElastic) {
            DyfStaggeredXSparse.redistribute(dist, dist);
            DyfStaggeredZSparse.redistribute(dist, dist);
            
            if(useFreeSurface != 1) {// if elastic && varGrid && useFreesurface==1, DybStaggeredXSparse nad DybStaggeredZSparse are not used
            DybStaggeredXSparse.redistribute(dist, dist);
            DybStaggeredZSparse.redistribute(dist, dist);
            }
            
            InterpolationStaggeredX.redistribute(dist, dist);
            InterpolationStaggeredZ.redistribute(dist, dist);
            InterpolationStaggeredXZ.redistribute(dist, dist);
        }

        InterpolationFull.redistribute(dist, dist);
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
    
    if ((useFreeSurface == 1) && (isElastic) && (!useVarGrid)){
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
        // + interpolation for Sxx,Syy,Szz,P
        numInterpMatrices += 1;
        if (isElastic) {
            // + Interpolation for Vx,Vz,Sxz
            numInterpMatrices += 3;
        }
    }


    return numInterpMatrices;
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
   if(isSetup){
        HOST_PRINT(dist->getCommunicatorPtr(), "", "Warning from estimateMemory FDTD3D.cpp : Existing Derivative Setup will be overwritten\n");
    }
    this->setup(config);   
    
return(estimateMemory(dist,modelCoordinates));
    
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::estimateMemory(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{

 SCAI_ASSERT(isSetup, "call setup function before estimateMemory");
return(this->getMemoryUsage(dist,modelCoordinates,getNumDMatrices(),getNumInterpMatrices()));
    
}

//! \brief Initializsation of the derivative matrices
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

    Dxf.scale(this->DT / DH);
    Dzf.scale(this->DT / DH);
    Dxb.scale(this->DT / DH);
    Dzb.scale(this->DT / DH);
    Dyf.scale(this->DT / DH);
    Dyb.scale(this->DT / DH);

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
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{

    SCAI_REGION("initializeMatrices")

    HOST_PRINT(comm, "", "Initialization of the matrices Dxf, Dyf, Dzf, Dxb, Dyb, Dzb\n");

    this->calcDxf(modelCoordinates, dist);
    this->calcDyf(modelCoordinates, dist);
    this->calcDzf(modelCoordinates, dist);
    this->calcDxb(modelCoordinates, dist);
    this->calcDyb(modelCoordinates, dist);
    this->calcDzb(modelCoordinates, dist);


    DxfSparse.setContextPtr(ctx);
    DxbSparse.setContextPtr(ctx);
    DyfSparse.setContextPtr(ctx);
    DybSparse.setContextPtr(ctx);
    DzfSparse.setContextPtr(ctx);
    DzbSparse.setContextPtr(ctx);

    DxfSparse.scale(this->DT);
    DxbSparse.scale(this->DT);
    DyfSparse.scale(this->DT);
    DybSparse.scale(this->DT);
    DzfSparse.scale(this->DT);
    DzbSparse.scale(this->DT);

    if ((isElastic) && (useVarGrid)) {
        this->calcDyfStaggeredX(modelCoordinates, dist);
        this->calcDybStaggeredX(modelCoordinates, dist);
        this->calcDyfStaggeredZ(modelCoordinates, dist);
        this->calcDybStaggeredZ(modelCoordinates, dist);

        DyfStaggeredXSparse.setContextPtr(ctx);
        DybStaggeredXSparse.setContextPtr(ctx);
        DyfStaggeredZSparse.setContextPtr(ctx);
        DybStaggeredZSparse.setContextPtr(ctx);

        DyfStaggeredXSparse *= this->DT;
        DybStaggeredXSparse *= this->DT;
        DyfStaggeredZSparse *= this->DT;
        DybStaggeredZSparse *= this->DT;
    }

    if (useVarGrid) {
        this->calcInterpolationFull(modelCoordinates, dist);
        if (isElastic) {
            this->calcInterpolationStaggeredX(modelCoordinates, dist);
            this->calcInterpolationStaggeredZ(modelCoordinates, dist);
            this->calcInterpolationStaggeredXZ(modelCoordinates, dist);
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
void KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::initializeFreeSurfaceMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
{


    this->calcDyfFreeSurface(modelCoordinates, dist);
    this->getDyfFreeSurface().setContextPtr(ctx);
    this->getDyfFreeSurface() *= this->DT;

    if (isElastic) {
        if (!useVarGrid) {
            this->calcDybFreeSurface(modelCoordinates, dist);
            this->getDybFreeSurface().setContextPtr(ctx);
            this->getDybFreeSurface() *= this->DT;

        } else {

            this->calcDybStaggeredXFreeSurface(modelCoordinates, dist);
            this->calcDybStaggeredZFreeSurface(modelCoordinates, dist);
            DybStaggeredXFreeSurface.setContextPtr(ctx);
            DybStaggeredZFreeSurface.setContextPtr(ctx);
            DybStaggeredXFreeSurface *= this->DT;
            DybStaggeredZFreeSurface *= this->DT;
            DybStaggeredXSparse.purge(); // DybStaggeredX won't be used with  varGrid+FreeSurface
            DybStaggeredZSparse.purge(); // DybStaggeredZ won't be used with  varGrid+FreeSurface
            DyfSparse.purge(); // DyfSparse won't be used with varGrid+FreeSurface
        }
    } else {
        //In acoustic modelling Dyf wont be used when free surface is used (see forwardsolver 2D/3D acoustic)
        this->getDyf().purge();
    }
    
}

//! \brief Getter method for combined derivative matrix
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::ForwardSolver::Derivatives::FDTD3D<ValueType>::getCombinedMatrix()
{
    auto temp = DxbSparse;
    temp+= DybSparse;
    temp+= DzbSparse;
    temp= transpose(temp);
    temp-=DxbSparse;
    temp-=DybSparse;
    temp-=DzbSparse;
//     auto temp2 = DxfSparse;
//     temp2 = transpose(temp);
//     temp-=temp2;
    return (temp);
}

template class KITGPI::ForwardSolver::Derivatives::FDTD3D<double>;
template class KITGPI::ForwardSolver::Derivatives::FDTD3D<float>;
