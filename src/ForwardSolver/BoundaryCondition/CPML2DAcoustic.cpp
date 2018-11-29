#include "CPML2DAcoustic.hpp"
using namespace scai;

//! \brief resetting the CPML memory variables
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::resetCPML()
{
    this->resetVector(psi_vxx);
    this->resetVector(psi_vyy);

    this->resetVector(psi_p_x);
    this->resetVector(psi_p_y);
}

//! \brief application of cpml on the derivation of vx in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::apply_vxx(scai::lama::Vector<ValueType> &vxx)
{
    this->applyCPML(vxx, psi_vxx, a_x, b_x, k_x);
}

//! \brief application of cpml on the derivation of vy in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::apply_vyy(scai::lama::Vector<ValueType> &vyy)
{
    this->applyCPML(vyy, psi_vyy, a_y, b_y, k_y);
}

//! \brief application of cpml on the derivation of p in x direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::apply_p_x(scai::lama::Vector<ValueType> &p_x)
{
    this->applyCPML(p_x, psi_p_x, a_x_half, b_x_half, k_x_half);
}

//! \brief application of cpml on the derivation of p in y direction
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::apply_p_y(scai::lama::Vector<ValueType> &p_y)
{
    this->applyCPML(p_y, psi_p_y, a_y_half, b_y_half, k_y_half);
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
 \param useFreeSurface Indicator which free surface is in use
 \param NPower degree of the damping profile
 \param KMaxCPML 
 \param CenterFrequencyCPML Center frequency inside the boundaries
 \param VMaxCPML Maximum p-wave velocity in the boundaries
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<ValueType>::init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, scai::IndexType useFreeSurface)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the PMl Coefficients...\n");
    
    active = true;

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

    this->initVector(psi_vxx, ctx, dist);
    this->initVector(psi_vyy, ctx, dist);

    this->initVector(psi_p_x, ctx, dist);
    this->initVector(psi_p_y, ctx, dist);

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

    Acquisition::Coordinates coordTransform(NX,NY,NZ);
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        coordinate = coordTransform.index2coordinate(read_localIndices_temp);
        gdist = coordTransform.edgeDistance(coordinate);

        if ((gdist.x < BoundaryWidth) || (gdist.y < BoundaryWidth)) {
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

template class KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML2DAcoustic<double>;
