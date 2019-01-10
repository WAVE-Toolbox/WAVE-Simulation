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
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Time sampling
 \param DH Grid spacing
 \param BoundaryWidth Width of damping boundary
 \param useFreeSurface Indicator which free surface is in use
 \param NPower degree of the damping profile
 \param KMaxCPML
 \param CenterFrequencyCPML Center frequency inside the boundaries
 \param VMaxCPML Maximum p-wave velocity in the boundaries
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, scai::IndexType useFreeSurface)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the PMl Coefficients...\n");

    active = true;

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

    /* Distributed vectors */
    k_x.setSameValue(dist, 1.0);
    k_x_half.setSameValue(dist, 1.0);
    k_y.setSameValue(dist, 1.0);
    k_y_half.setSameValue(dist, 1.0);
    k_z.setSameValue(dist, 1.0);
    k_z_half.setSameValue(dist, 1.0);

    lama::DenseVector<ValueType> k_x_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_x_half_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_y_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_y_half_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_z_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_z_half_temp(dist, 1.0, ctx);

    lama::DenseVector<ValueType> a_x_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_x_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_y_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_y_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_z_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_z_half_temp(dist, 0.0, ctx);

    lama::DenseVector<ValueType> b_x_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_x_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_y_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_y_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_z_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_z_half_temp(dist, 0.0, ctx);

    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = modelCoordinates.index2coordinate(read_localIndices_temp);
        gdist = modelCoordinates.edgeDistance(coordinate);

        if (gdist.min() < BoundaryWidth) {
            if (gdist.x < BoundaryWidth) {
                this->SetCoeffCPML(a_x_temp, b_x_temp, k_x_temp, a_x_half_temp, b_x_half_temp, k_x_half_temp, coordinate.x, gdist.x, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
            }
            if (gdist.y < BoundaryWidth) {
                this->SetCoeffCPML(a_y_temp, b_y_temp, k_y_temp, a_y_half_temp, b_y_half_temp, k_y_half_temp, coordinate.y, gdist.y, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
                if (useFreeSurface > 0) {
                    if (coordinate.y < BoundaryWidth) {
                        this->ResetCoeffFreeSurface(a_y_temp, b_y_temp, k_y_temp, a_y_half_temp, b_y_half_temp, k_y_half_temp, i);
                    }
                }
            }
            if (gdist.z < BoundaryWidth) {
                this->SetCoeffCPML(a_z_temp, b_z_temp, k_z_temp, a_z_half_temp, b_z_half_temp, k_z_half_temp, coordinate.z, gdist.z, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
            }
        }
    }

    k_x = k_x_temp;
    k_x_half = k_x_half_temp;
    k_y = k_y_temp;
    k_y_half = k_y_half_temp;
    k_z = k_z_temp;
    k_z_half = k_z_half_temp;

    a_x = a_x_temp;
    a_x_half = a_x_half_temp;
    a_y = a_y_temp;
    a_y_half = a_y_half_temp;
    a_z = a_z_temp;
    a_z_half = a_z_half_temp;

    b_x = b_x_temp;
    b_x_half = b_x_half_temp;
    b_y = b_y_temp;
    b_y_half = b_y_half_temp;
    b_z = b_z_temp;
    b_z_half = b_z_half_temp;

    //     /* Release all read and write access */
    read_localIndices.release();

    HOST_PRINT(comm, "", "Finished with initialization of the CPML coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML3DAcoustic<double>;
