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
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vxx(scai::lama::Vector<ValueType> &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vyy(scai::lama::Vector<ValueType> &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of vz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_vzz(scai::lama::Vector<ValueType> &vzz)
{
    this->applyCPML(vzz, psi_vzz, a_z, b_z, k_z);
}

//! \brief application of cpml on the derivation of p in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_x(scai::lama::Vector<ValueType> &p_x)
{
    this->applyCPML(p_x, psi_p_x, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of p in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_y(scai::lama::Vector<ValueType> &p_y)
{
    this->applyCPML(p_y, psi_p_y, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of p in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::apply_p_z(scai::lama::Vector<ValueType> &p_z)
{
    this->applyCPML(p_z, psi_p_z, a_z_half, b_z_half, k_z_half);
}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param DT Time sampling
 \param DH Grid spacing
 \param BoundaryWidth Width of damping boundary
 \param useFreeSurface Bool if free surface is in use
 \param NPower degree of the damping profile
 \param KMaxCPML
 \param CenterFrequencyCPML Center frequency inside the boundaries
 \param VMaxCPML Maximum p-wave velocity in the boundaries
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, scai::IndexType useFreeSurface)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the PMl Coefficients...\n");

    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);

    IndexType numLocalIndices = localIndices.size(); // Number of local indices

    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp;                             // Temporary storage, so we do not have to access the array

    // 	update_PmlTemp.allocate ( dist );
    // 	update_PmlTemp=0;

    /* Distributed vectors */

    /* Distributed vectors */
    this->initVector(psi_vxx, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);
    this->initVector(psi_vzz, ctx, dist);

    this->initVector(psi_p_x, ctx, dist);
    this->initVector(psi_p_y, ctx, dist);
    this->initVector(psi_p_z, ctx, dist);

    this->initVector(k_x, ctx, dist);
    this->initVector(k_y, ctx, dist);
    this->initVector(k_z, ctx, dist);
    this->initVector(b_x, ctx, dist);
    this->initVector(b_y, ctx, dist);
    this->initVector(b_z, ctx, dist);
    this->initVector(a_x, ctx, dist);
    this->initVector(a_y, ctx, dist);
    this->initVector(a_z, ctx, dist);

    this->initVector(k_x_half, ctx, dist);
    this->initVector(k_y_half, ctx, dist);
    this->initVector(k_z_half, ctx, dist);
    this->initVector(b_x_half, ctx, dist);
    this->initVector(b_y_half, ctx, dist);
    this->initVector(b_z_half, ctx, dist);
    this->initVector(a_x_half, ctx, dist);
    this->initVector(a_y_half, ctx, dist);
    this->initVector(a_z_half, ctx, dist);

    k_x = 1.0;
    k_y = 1.0;
    k_z = 1.0;
    k_x_half = 1.0;
    k_y_half = 1.0;
    k_z_half = 1.0;

    Acquisition::Coordinates coordTransform;
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ);
        gdist = coordTransform.edgeDistance(coordinate, NX, NY, NZ);

        if (gdist.min() < BoundaryWidth) {
            if (gdist.x < BoundaryWidth) {
                this->SetCoeffCPML(a_x, b_x, k_x, a_x_half, b_x_half, k_x_half, coordinate.x, gdist.x, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
            }
            if (gdist.y < BoundaryWidth) {
                this->SetCoeffCPML(a_y, b_y, k_y, a_y_half, b_y_half, k_y_half, coordinate.y, gdist.y, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
                if (useFreeSurface > 0) {
                    if (coordinate.y < BoundaryWidth) {
                        this->ResetCoeffFreeSurface(a_y, b_y, k_y, a_y_half, b_y_half, k_y_half, i);
                    }
                }
            }
            if (gdist.z < BoundaryWidth) {
                this->SetCoeffCPML(a_z, b_z, k_z, a_z_half, b_z_half, k_z_half, coordinate.z, gdist.z, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
            }
        }
    }
    //

    //     /* Release all read and write access */
    read_localIndices.release();

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished with initialization of the CPML coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<double>;
