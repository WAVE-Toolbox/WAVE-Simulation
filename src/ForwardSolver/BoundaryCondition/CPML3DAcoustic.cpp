#include "CPML3DAcoustic.hpp"
using namespace scai;

//! \brief resetting the CPML memory variables
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::resetCPML()
{
    this->resetVector(psi_vxx);
    this->resetVector(psi_vyy);
    this->resetVector(psi_vzz);

    this->resetVector(psi_p_x);
    this->resetVector(psi_p_y);
    this->resetVector(psi_p_z);
}

//! \brief application of cpml on the derivation of vx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vxx(scai::lama::DenseVector<ValueType> &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vyy(scai::lama::DenseVector<ValueType> &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y);
}

//! \brief application of cpml on the derivation of vz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vzz(scai::lama::DenseVector<ValueType> &vzz)
{
    this->applyCPML(vzz, psi_vzz, a_z, b_z);
}

//! \brief application of cpml on the derivation of p in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_x(scai::lama::DenseVector<ValueType> &p_x)
{
    this->applyCPML(p_x, psi_p_x, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of p in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_y(scai::lama::DenseVector<ValueType> &p_y)
{
    this->applyCPML(p_y, psi_p_y, a_y_half, b_y_half);
}

//! \brief application of cpml on the derivation of p in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_z(scai::lama::DenseVector<ValueType> &p_z)
{
    this->applyCPML(p_z, psi_p_z, a_z_half, b_z_half);
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
ValueType KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
        IndexType zDist = gdist.z / modelCoordinates.getDHFactor(layer);

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

        if (zDist < width) {
            counter++;
        }
    }
    IndexType sum=dist->getCommunicator().sum(counter);
    IndexType numVectorsPerDim = 6;
    return(sum * sizeof(ValueType) * numVectorsPerDim / (1024 * 1024));
}

//! \brief Initialization of the absorbing coefficient matrix
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
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the PMl Coefficients...\n");

    active = true;

    /* Distributed vectors */
    this->initVector(psi_vxx, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);
    this->initVector(psi_vzz, ctx, dist);

    this->initVector(psi_p_x, ctx, dist);
    this->initVector(psi_p_y, ctx, dist);
    this->initVector(psi_p_z, ctx, dist);

    this->initVector(a_x, ctx, dist);
    this->initVector(b_x, ctx, dist);
    this->initVector(a_x_half, ctx, dist);
    this->initVector(b_x_half, ctx, dist);

    this->initVector(a_y, ctx, dist);
    this->initVector(b_y, ctx, dist);
    this->initVector(a_y_half, ctx, dist);
    this->initVector(b_y_half, ctx, dist);

    this->initVector(a_z, ctx, dist);
    this->initVector(b_z, ctx, dist);
    this->initVector(a_z_half, ctx, dist);
    this->initVector(b_z_half, ctx, dist);

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
    lama::VectorAssembly<ValueType> a_zAssembly, b_zAssembly, a_z_halfAssembly, b_z_halfAssembly;

    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        coordinate = modelCoordinates.index2coordinate(ownedIndex);
        IndexType layer = modelCoordinates.getLayer(coordinate);

        gdist = modelCoordinates.edgeDistance(coordinate);

        IndexType width = std::ceil((float)BoundaryWidth / modelCoordinates.getDHFactor(layer));
        IndexType xDist = gdist.x / modelCoordinates.getDHFactor(layer);
        IndexType yDist = gdist.y / modelCoordinates.getDHFactor(layer);
        IndexType zDist = gdist.z / modelCoordinates.getDHFactor(layer);

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

        if (zDist < width) {
            IndexType zCoord = coordinate.z / modelCoordinates.getDHFactor(layer);
            if (zCoord < width) {
                a_zAssembly.push(ownedIndex, a.at(layer)[zDist]);
                b_zAssembly.push(ownedIndex, b.at(layer)[zDist]);
                a_z_halfAssembly.push(ownedIndex, a_half.at(layer)[zDist]);
                b_z_halfAssembly.push(ownedIndex, b_half.at(layer)[zDist]);
            } else {
                a_zAssembly.push(ownedIndex, a_half.at(layer)[zDist]);
                b_zAssembly.push(ownedIndex, b_half.at(layer)[zDist]);
                a_z_halfAssembly.push(ownedIndex, a.at(layer)[zDist]);
                b_z_halfAssembly.push(ownedIndex, b.at(layer)[zDist]);
            }
        }
    }

    a_x.fillFromAssembly(a_xAssembly);
    b_x.fillFromAssembly(b_xAssembly);
    a_x_half.fillFromAssembly(a_x_halfAssembly);
    b_x_half.fillFromAssembly(b_x_halfAssembly);

    a_y.fillFromAssembly(a_yAssembly);
    b_y.fillFromAssembly(b_yAssembly);
    a_y_half.fillFromAssembly(a_y_halfAssembly);
    b_y_half.fillFromAssembly(b_y_halfAssembly);

    a_z.fillFromAssembly(a_zAssembly);
    b_z.fillFromAssembly(b_zAssembly);
    a_z_half.fillFromAssembly(a_z_halfAssembly);
    b_z_half.fillFromAssembly(b_z_halfAssembly);

    HOST_PRINT(comm, "", "Finished with initialization of the CPML coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<double>;
