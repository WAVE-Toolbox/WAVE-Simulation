#include "CPML3D.hpp"
using namespace scai;

//! \brief resetting the CPML memory variables
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::reset()
{
    this->resetVector(psi_vxx);
    this->resetVector(psi_vyx);
    this->resetVector(psi_vzx);
    this->resetVector(psi_vxy);
    this->resetVector(psi_vyy);
    this->resetVector(psi_vzy);
    this->resetVector(psi_vxz);
    this->resetVector(psi_vyz);
    this->resetVector(psi_vzz);

    this->resetVector(psi_sxx_x);
    this->resetVector(psi_sxy_x);
    this->resetVector(psi_sxz_x);
    this->resetVector(psi_sxy_y);
    this->resetVector(psi_syy_y);
    this->resetVector(psi_syz_y);
    this->resetVector(psi_sxz_z);
    this->resetVector(psi_syz_z);
    this->resetVector(psi_szz_z);
}

//! \brief application of cpml on the derivation of sxx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_sxx_x(scai::lama::Vector &sxx_x)
{
    this->applyCPML(sxx_x, psi_sxx_x, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of sxy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_sxy_x(scai::lama::Vector &sxy_x)
{
    this->applyCPML(sxy_x, psi_sxy_x, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of sxz in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_sxz_x(scai::lama::Vector &sxz_x)
{
    this->applyCPML(sxz_x, psi_sxz_x, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of sxy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_sxy_y(scai::lama::Vector &sxy_y)
{
    this->applyCPML(sxy_y, psi_sxy_y, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of syy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_syy_y(scai::lama::Vector &syy_y)
{
    this->applyCPML(syy_y, psi_syy_y, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of syz in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_syz_y(scai::lama::Vector &syz_y)
{
    this->applyCPML(syz_y, psi_syz_y, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of sxz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_sxz_z(scai::lama::Vector &sxz_z)
{
    this->applyCPML(sxz_z, psi_sxz_z, a_z, b_z, k_z);
}

//! \brief application of cpml on the derivation of syz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_syz_z(scai::lama::Vector &syz_z)
{
    this->applyCPML(syz_z, psi_syz_z, a_z, b_z, k_z);
}

//! \brief application of cpml on the derivation of szz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_szz_z(scai::lama::Vector &szz_z)
{
    this->applyCPML(szz_z, psi_szz_z, a_z_half, b_z_half, k_z_half);
}

//! \brief application of cpml on the derivation of vx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vxx(scai::lama::Vector &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of vy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vyx(scai::lama::Vector &vyx)
{
    this->applyCPML(vyx, psi_vyx, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of vz in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vzx(scai::lama::Vector &vzx)
{
    this->applyCPML(vzx, psi_vzx, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of vx in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vxy(scai::lama::Vector &vxy)
{
    this->applyCPML(vxy, psi_vxy, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vyy(scai::lama::Vector &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of vz in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vzy(scai::lama::Vector &vzy)
{
    this->applyCPML(vzy, psi_vzy, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of vx in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vxz(scai::lama::Vector &vxz)
{
    this->applyCPML(vxz, psi_vxz, a_z_half, b_z_half, k_z_half);
}

//! \brief application of cpml on the derivation of vy in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vyz(scai::lama::Vector &vyz)
{
    this->applyCPML(vyz, psi_vyz, a_z_half, b_z_half, k_z_half);
}

//! \brief application of cpml on the derivation of vz in z direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::apply_vzz(scai::lama::Vector &vzz)
{
    this->applyCPML(vzz, psi_vzz, a_z, b_z, k_z);
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
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface)
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
    this->initVector(psi_vyx, ctx, dist);
    this->initVector(psi_vzx, ctx, dist);
    this->initVector(psi_vxy, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);
    this->initVector(psi_vzy, ctx, dist);
    this->initVector(psi_vxz, ctx, dist);
    this->initVector(psi_vyz, ctx, dist);
    this->initVector(psi_vzz, ctx, dist);

    this->initVector(psi_sxx_x, ctx, dist);
    this->initVector(psi_sxy_x, ctx, dist);
    this->initVector(psi_sxz_x, ctx, dist);
    this->initVector(psi_sxy_y, ctx, dist);
    this->initVector(psi_syy_y, ctx, dist);
    this->initVector(psi_syz_y, ctx, dist);
    this->initVector(psi_sxz_z, ctx, dist);
    this->initVector(psi_syz_z, ctx, dist);
    this->initVector(psi_szz_z, ctx, dist);

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

    Acquisition::Coordinates<ValueType> coordTransform;
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
                if (useFreeSurface) {
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

template class KITGPI::ForwardSolver::BoundaryCondition::CPML3D<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML3D<double>;