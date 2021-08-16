#include "CPMLEM3D.hpp"
using namespace scai;

//! \brief resetting the CPMLEM memory variables
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::resetCPML()
{
    this->resetVector(psi_hyx);
    this->resetVector(psi_hzx);
    this->resetVector(psi_hxy);
    this->resetVector(psi_hzy);
    this->resetVector(psi_hxz);
    this->resetVector(psi_hyz);

    this->resetVector(psi_ezx);
    this->resetVector(psi_eyx);
    this->resetVector(psi_ezy);
    this->resetVector(psi_exy);
    this->resetVector(psi_eyz);
    this->resetVector(psi_exz);
}

//! \brief application of cpml on the derivation of sxy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_ezx(scai::lama::DenseVector<ValueType> &ezx)
{
    this->applyCPMLEM(ezx, psi_ezx, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of sxz in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_eyx(scai::lama::DenseVector<ValueType> &eyx)
{
    this->applyCPMLEM(eyx, psi_eyx, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of sxy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_ezy(scai::lama::DenseVector<ValueType> &ezy)
{
    this->applyCPMLEM(ezy, psi_ezy, a_y_half, b_y_half);
}

//! \brief application of cpml on the derivation of syz in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_exy(scai::lama::DenseVector<ValueType> &exy)
{
    this->applyCPMLEM(exy, psi_exy, a_y_half, b_y_half);
}

//! \brief application of cpml on the derivation of sxz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_eyz(scai::lama::DenseVector<ValueType> &eyz)
{
    this->applyCPMLEM(eyz, psi_eyz, a_z_half, b_z_half);
}

//! \brief application of cpml on the derivation of syz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_exz(scai::lama::DenseVector<ValueType> &exz)
{
    this->applyCPMLEM(exz, psi_exz, a_z_half, b_z_half);
}

//! \brief application of cpml on the derivation of hy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hyx(scai::lama::DenseVector<ValueType> &hyx)
{
    this->applyCPMLEM(hyx, psi_hyx, a_x_half, b_x_half);
}

//! \brief application of cpml on the derivation of hz in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hzx(scai::lama::DenseVector<ValueType> &hzx)
{
    this->applyCPMLEM(hzx, psi_hzx, a_x, b_x);
}

//! \brief application of cpml on the derivation of hx in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hxy(scai::lama::DenseVector<ValueType> &hxy)
{
    this->applyCPMLEM(hxy, psi_hxy, a_y, b_y);
}

//! \brief application of cpml on the derivation of hz in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hzy(scai::lama::DenseVector<ValueType> &hzy)
{
    this->applyCPMLEM(hzy, psi_hzy, a_y, b_y);
}

//! \brief application of cpml on the derivation of hx in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hxz(scai::lama::DenseVector<ValueType> &hxz)
{
    this->applyCPMLEM(hxz, psi_hxz, a_z, b_z);
}

//! \brief application of cpml on the derivation of hy in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::apply_hyz(scai::lama::DenseVector<ValueType> &hyz)
{
    this->applyCPMLEM(hyz, psi_hyz, a_z, b_z);
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
ValueType KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::estimateMemory(IndexType BoundaryWidth, IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
    IndexType numVectorsPerDim = 10;
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
 \param useFreeSurface Bool if free surface is in use
 \param NPower degree of the damping profile
 \param CenterFrequencyCPML Center frequency inside the boundaries
 \param VMaxCPML Maximum p-wave velocity in the boundaries
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<ValueType>::init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the PMl Coefficients...\n");

    active = true;

    /* Distributed vectors */
    this->initVector(psi_hyx, ctx, dist);
    this->initVector(psi_hzx, ctx, dist);
    this->initVector(psi_hxy, ctx, dist);
    this->initVector(psi_hzy, ctx, dist);
    this->initVector(psi_hxz, ctx, dist);
    this->initVector(psi_hyz, ctx, dist);

    this->initVector(psi_ezx, ctx, dist);
    this->initVector(psi_eyx, ctx, dist);
    this->initVector(psi_ezy, ctx, dist);
    this->initVector(psi_exy, ctx, dist);
    this->initVector(psi_eyz, ctx, dist);
    this->initVector(psi_exz, ctx, dist);

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

        this->calcCoeffCPMLEM(aLayer, bLayer, NPower, CenterFrequencyCPML, VMaxCPML, DT, modelCoordinates.getDH(layer));
        this->calcCoeffCPMLEM(a_halfLayer, b_halfLayer, NPower, CenterFrequencyCPML, VMaxCPML, DT, modelCoordinates.getDH(layer), 1);

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

    HOST_PRINT(comm, "", "Finished with initialization of the CPMLEM coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPMLEM3D<double>;
