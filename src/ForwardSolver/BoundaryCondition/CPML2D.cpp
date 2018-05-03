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
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxx_x(scai::lama::Vector<ValueType> &sxx_x)
{
    this->applyCPML(sxx_x, psi_sxx_x, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of sxy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxy_x(scai::lama::Vector<ValueType> &sxy_x)
{
    this->applyCPML(sxy_x, psi_sxy_x, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of sxy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_sxy_y(scai::lama::Vector<ValueType> &sxy_y)
{
    this->applyCPML(sxy_y, psi_sxy_y, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of syy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_syy_y(scai::lama::Vector<ValueType> &syy_y)
{
    this->applyCPML(syy_y, psi_syy_y, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of vx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vxx(scai::lama::Vector<ValueType> &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of vy in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vyx(scai::lama::Vector<ValueType> &vyx)
{
    this->applyCPML(vyx, psi_vyx, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of vx in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vxy(scai::lama::Vector<ValueType> &vxy)
{
    this->applyCPML(vxy, psi_vxy, a_y_half, b_y_half, k_y_half);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::apply_vyy(scai::lama::Vector<ValueType> &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y, k_y);
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
void KITGPI::ForwardSolver::BoundaryCondition::CPML2D<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface)
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
    this->initVector(psi_vxy, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);

    this->initVector(psi_sxx_x, ctx, dist);
    this->initVector(psi_sxy_x, ctx, dist);
    this->initVector(psi_sxy_y, ctx, dist);
    this->initVector(psi_syy_y, ctx, dist);

    std::cout << psi_syy_y << std::endl;

    /* Distributed vectors */
    k_x.setSameValue(dist, 1.0);
    k_y.setSameValue(dist, 1.0);
    k_x_half.setSameValue(dist, 1.0);
    k_y_half.setSameValue(dist, 1.0);

    lama::DenseVector<ValueType> k_x_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_y_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_x_half_temp(dist, 1.0, ctx);
    lama::DenseVector<ValueType> k_y_half_temp(dist, 1.0, ctx);

    lama::DenseVector<ValueType> a_x_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_y_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_x_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> a_y_half_temp(dist, 0.0, ctx);

    lama::DenseVector<ValueType> b_x_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_y_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_x_half_temp(dist, 0.0, ctx);
    lama::DenseVector<ValueType> b_y_half_temp(dist, 0.0, ctx);

    Acquisition::Coordinates coordTransform;
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ);
        gdist = coordTransform.edgeDistance(coordinate, NX, NY, NZ);

        if ((gdist.x < BoundaryWidth) || (gdist.y < BoundaryWidth)) {
            if (gdist.x < BoundaryWidth) {
                this->SetCoeffCPML(a_x_temp, b_x_temp, k_x_temp, a_x_half_temp, b_x_half_temp, k_x_half_temp, coordinate.x, gdist.x, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
            }
            if (gdist.y < BoundaryWidth) {
                this->SetCoeffCPML(a_y_temp, b_y_temp, k_y_temp, a_y_half_temp, b_y_half_temp, k_y_half_temp, coordinate.y, gdist.y, BoundaryWidth, NPower, KMaxCPML, CenterFrequencyCPML, VMaxCPML, i, DT, DH);
                if (useFreeSurface) {
                    if (coordinate.y < BoundaryWidth) {
                        this->ResetCoeffFreeSurface(a_y_temp, b_y_temp, k_y_temp, a_y_half_temp, b_y_half_temp, k_y_half_temp, i);
                    }
                }
            }
        }
    }
    //
    k_x = k_x_temp;
    k_y = k_y_temp;
    k_x_half = k_x_half_temp;
    k_y_half = k_y_half_temp;
    a_x = a_x_temp;
    a_y = a_y_temp;
    a_x_half = a_x_half_temp;
    a_y_half = a_y_half_temp;
    b_x = b_x_temp;
    b_y = b_y_temp;
    b_x_half = b_x_half_temp;
    b_y_half = b_y_half_temp;

    //     /* Release all read and write access */
    read_localIndices.release();

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished with initialization of the CPML coefficients!\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML2D<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML2D<float>;
