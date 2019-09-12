#include "CPML2D.hpp"
using namespace scai;

//! \brief resetting the CPML memory variables
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::resetCPML()
{
    this->resetVector(psi_vxx);
    this->resetVector(psi_vyx);
    this->resetVector(psi_vxy);
    this->resetVector(psi_vyy);

    this->resetVector(psi_sxx_x);
    this->resetVector(psi_sxy_x);
    this->resetVector(psi_sxy_y);
    this->resetVector(psi_syy_y);
}

//! \brief application of cpml on the derivation of sxx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxx_x(scai::lama::DenseVector<ValueType> &sxx_x)
{
    this->applyCPML(sxx_x, psi_sxx_x, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of sxy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxy_x(scai::lama::DenseVector<ValueType> &sxy_x)
{
    this->applyCPML(sxy_x, psi_sxy_x, a_x, b_x);
}

//! \brief application of cpml on the derivation of sxy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxy_y(scai::lama::DenseVector<ValueType> &sxy_y)
{
    this->applyCPML(sxy_y, psi_sxy_y, a_y, b_y);
}

//! \brief application of cpml on the derivation of syy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_syy_y(scai::lama::DenseVector<ValueType> &syy_y)
{
    this->applyCPML(syy_y, psi_syy_y, a_y_half, b_y_half);
}

//! \brief application of cpml on the derivation of vx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vxx(scai::lama::DenseVector<ValueType> &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x);
}

//! \brief application of cpml on the derivation of vy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vyx(scai::lama::DenseVector<ValueType> &vyx)
{
    this->applyCPML(vyx, psi_vyx, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of vx in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vxy(scai::lama::DenseVector<ValueType> &vxy)
{
    this->applyCPML(vxy, psi_vxy, a_y_half, b_y_half);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vyy(scai::lama::DenseVector<ValueType> &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y);
}

//! \brief estimate memory for the absorbing boundary frame
/*!
 *
 \param BoundaryWidth Width of damping boundary
 \param useFreeSurface Indicator which free surface is in use
  \param dist Distribution of the wavefield
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{

    IndexType counter = 0;

    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        IndexType layer = modelCoordinates.getLayer(coordinate);

        Acquisition::coordinate3D gdist = modelCoordinates.edgeDistance(coordinate);

        IndexType width = std::ceil((float)BoundaryWidth / modelCoordinates.getDHFactor(layer));
        IndexType xDist = gdist.x / modelCoordinates.getDHFactor(layer);
        IndexType yDist = gdist.y / modelCoordinates.getDHFactor(layer);

        if (xDist < width) {
            counter++;
        }

        if (yDist < width) {
            IndexType yCoord = coordinate.y / modelCoordinates.getDHFactor(layer);
            if (yCoord < width) {
                if (useFreeSurface == 0) {
                    counter++;
                }
            } else {
                counter++;
            }
        }
    }

     IndexType sum=dist->getCommunicator().sum(counter);
    IndexType numVectorsPerDim = 8;
    return(sum * sizeof(ValueType) * numVectorsPerDim / (1024 * 1024));
}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Time sampling
 \param BoundaryWidth Width of damping boundary
 \param useFreeSurface Indicator which free surface is in use
 \param NPower degree of the damping profile
 \param CenterFrequencyCPML Center frequency inside the boundaries
 \param VMaxCPML Maximum p-wave velocity in the boundaries
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the PMl Coefficients...\n");

    active = true;

    /* Distributed vectors */
    this->initVector(psi_vxx, ctx, dist);
    this->initVector(psi_vyx, ctx, dist);
    this->initVector(psi_vxy, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);

    this->initVector(psi_sxx_x, ctx, dist);
    this->initVector(psi_sxy_x, ctx, dist);
    this->initVector(psi_sxy_y, ctx, dist);
    this->initVector(psi_syy_y, ctx, dist);

    this->initVector(a_x, ctx, dist);
    this->initVector(b_x, ctx, dist);
    this->initVector(a_x_half, ctx, dist);
    this->initVector(b_x_half, ctx, dist);

    this->initVector(a_y, ctx, dist);
    this->initVector(b_y, ctx, dist);
    this->initVector(a_y_half, ctx, dist);
    this->initVector(b_y_half, ctx, dist);

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

    std::vector<std::vector<ValueType>> a, b, a_half, b_half;

    for (IndexType layer = 0; layer < modelCoordinates.getNumLayers(); layer++) {

        IndexType width = std::ceil((float)BoundaryWidth / modelCoordinates.getDHFactor(layer));
        std::vector<ValueType> aLayer(width, 0.0), bLayer(width, 0.0), a_halfLayer(width, 0.0), b_halfLayer(width, 0.0);

        this->calcCoeffCPML(aLayer, bLayer, NPower, CenterFrequencyCPML, VMaxCPML, DT, modelCoordinates.getDH(layer));
        this->calcCoeffCPML(a_halfLayer, b_halfLayer, NPower, CenterFrequencyCPML, VMaxCPML, DT, modelCoordinates.getDH(layer), 1);

        a.push_back(aLayer);
        b.push_back(bLayer);
        a_half.push_back(a_halfLayer);
        b_half.push_back(b_halfLayer);
    }

    lama::VectorAssembly<ValueType> a_xAssembly, b_xAssembly, a_x_halfAssembly, b_x_halfAssembly;
    lama::VectorAssembly<ValueType> a_yAssembly, b_yAssembly, a_y_halfAssembly, b_y_halfAssembly;

    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        coordinate = modelCoordinates.index2coordinate(ownedIndex);
        IndexType layer = modelCoordinates.getLayer(coordinate);

        gdist = modelCoordinates.edgeDistance(coordinate);

        IndexType width = std::ceil((float)BoundaryWidth / modelCoordinates.getDHFactor(layer));
        IndexType xDist = gdist.x / modelCoordinates.getDHFactor(layer);
        IndexType yDist = gdist.y / modelCoordinates.getDHFactor(layer);

        if (xDist < width) {
            IndexType xCoord = coordinate.x / modelCoordinates.getDHFactor(layer);
            if (xCoord < width) {
                a_xAssembly.push(ownedIndex, a.at(layer)[xDist]);
                b_xAssembly.push(ownedIndex, b.at(layer)[xDist]);
                a_x_halfAssembly.push(ownedIndex, a_half.at(layer)[xDist]);
                b_x_halfAssembly.push(ownedIndex, b_half.at(layer)[xDist]);
            } else {
                a_xAssembly.push(ownedIndex, a_half.at(layer)[xDist]);
                b_xAssembly.push(ownedIndex, b_half.at(layer)[xDist]);
                a_x_halfAssembly.push(ownedIndex, a.at(layer)[xDist]);
                b_x_halfAssembly.push(ownedIndex, b.at(layer)[xDist]);
            }
        }

        if (yDist < width) {
            IndexType yCoord = coordinate.y / modelCoordinates.getDHFactor(layer);
            if (yCoord < width) {
                if (useFreeSurface == 0) {
                    a_yAssembly.push(ownedIndex, a.at(layer)[yDist]);
                    b_yAssembly.push(ownedIndex, b.at(layer)[yDist]);
                    a_y_halfAssembly.push(ownedIndex, a_half.at(layer)[yDist]);
                    b_y_halfAssembly.push(ownedIndex, b_half.at(layer)[yDist]);
                }
            } else {
                a_yAssembly.push(ownedIndex, a_half.at(layer)[yDist]);
                b_yAssembly.push(ownedIndex, b_half.at(layer)[yDist]);
                a_y_halfAssembly.push(ownedIndex, a.at(layer)[yDist]);
                b_y_halfAssembly.push(ownedIndex, b.at(layer)[yDist]);
            }
        }
    }

    a_x.fillFromAssembly(a_xAssembly);
    a_y.fillFromAssembly(a_yAssembly);
    a_x_half.fillFromAssembly(a_x_halfAssembly);
    a_y_half.fillFromAssembly(a_y_halfAssembly);
    b_x.fillFromAssembly(b_xAssembly);
    b_y.fillFromAssembly(b_yAssembly);
    b_x_half.fillFromAssembly(b_x_halfAssembly);
    b_y_half.fillFromAssembly(b_y_halfAssembly);

    HOST_PRINT(comm, "", "Finished with initialization of the CPML coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML2D<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML2D<float>;
