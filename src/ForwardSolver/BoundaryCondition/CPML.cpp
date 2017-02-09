#include "CPML.hpp"
using namespace scai;

/*! \brief Reset a single Vector to zero.
*/
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::resetVector(scai::lama::DenseVector<ValueType> &vector)
{
    vector = 0.0;
}

/*! \brief Intitialisation of a single vector.
*
* This method will set the context, allocate the the wavefield and set the field to zero.
*/
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::initVector(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);

    resetVector(vector);
}

/*! \brief set CPML coefficients
 * 
 * method to set cpml coefficients for a given gridpoint
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::SetCoeffCPML(scai::lama::DenseVector<ValueType> &a, scai::lama::DenseVector<ValueType> &b, scai::lama::DenseVector<ValueType> &kInv, scai::lama::DenseVector<ValueType> &a_half, scai::lama::DenseVector<ValueType> &b_half, scai::lama::DenseVector<ValueType> &kInv_half, IndexType coord, IndexType gdist, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, IndexType i, ValueType DT, ValueType DH)
{

    ValueType RCoef = 0.0008;

    ValueType alpha_max_Pml = 2.0 * M_PI * (CenterFrequencyCPML / 2.0);

    ValueType d0 = -(NPower + 1) * VMaxCPML * log(RCoef) / (2.0 * BoundaryWidth * DH);

    ValueType PositionNorm = 0.0;

    ValueType k_temp = 0.0;
    ValueType b_temp = 0.0;
    ValueType a_temp = 0.0;

    ValueType alpha_prime = 0.0;
    ValueType d = 0.0;

    utilskernel::LArray<ValueType> *k_LA = &kInv.getLocalValues();
    hmemo::WriteAccess<ValueType> write_k(*k_LA);

    utilskernel::LArray<ValueType> *b_LA = &b.getLocalValues();
    hmemo::WriteAccess<ValueType> write_b(*b_LA);

    utilskernel::LArray<ValueType> *a_LA = &a.getLocalValues();
    hmemo::WriteAccess<ValueType> write_a(*a_LA);

    utilskernel::LArray<ValueType> *k_half_LA = &kInv_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_k_half(*k_half_LA);

    utilskernel::LArray<ValueType> *b_half_LA = &b_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_b_half(*b_half_LA);

    utilskernel::LArray<ValueType> *a_half_LA = &a_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_a_half(*a_half_LA);

    /* left boundary */
    if (coord < BoundaryWidth) {
        PositionNorm = (ValueType)(BoundaryWidth - gdist) / BoundaryWidth;
        d = d0 * pow(PositionNorm, NPower);
        k_temp = 1.0 + (KMaxCPML - 1.0) * pow(PositionNorm, NPower);
        alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
        b_temp = exp(-(d / k_temp + alpha_prime) * DT);
        /* avoid division by zero outside the PML */
        if (abs(d) > 1.0e-6) {
            a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));
        } else
            a_temp = 0.0;
        write_a[i] = a_temp;
        write_b[i] = b_temp;
        write_k[i] = 1 / k_temp;

        /* right boundary
		 *starts with half point -> first point shiftet +1 to th right*/
    } else if (gdist < BoundaryWidth - 1) {
        PositionNorm = (ValueType)(BoundaryWidth - gdist - 1) / BoundaryWidth;
        d = d0 * pow(PositionNorm, NPower);
        k_temp = 1.0 + (KMaxCPML - 1.0) * pow(PositionNorm, NPower);
        alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
        b_temp = exp(-(d / k_temp + alpha_prime) * DT);
        /* avoid division by zero outside the PML */
        if (abs(d) > 1.0e-6) {
            a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));
        } else
            a_temp = 0.0;
        write_a[i] = a_temp;
        write_b[i] = b_temp;
        write_k[i] = 1 / k_temp;
    }

    /* half points */
    PositionNorm = (ValueType)(BoundaryWidth - gdist - 0.5) / BoundaryWidth;
    d = d0 * pow(PositionNorm, NPower);
    k_temp = 1.0 + (KMaxCPML - 1.0) * pow(PositionNorm, NPower);
    alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
    b_temp = exp(-(d / k_temp + alpha_prime) * DT);
    /* avoid division by zero outside the PML */
    if (abs(d) > 1.0e-6) {
        a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));
    } else
        a_temp = 0.0;
    write_a_half[i] = a_temp;
    write_b_half[i] = b_temp;
    write_k_half[i] = 1 / k_temp;

    write_k.release();
    write_a.release();
    write_b.release();
    write_k_half.release();
    write_a_half.release();
    write_b_half.release();
}

/*! \brief Reset CPML coefficients for free surface
 * 
 * method to set cpml coefficients for a given gridpoint
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::ResetCoeffFreeSurface(lama::DenseVector<ValueType> &a, lama::DenseVector<ValueType> &b, lama::DenseVector<ValueType> &kInv,
                                                                                      lama::DenseVector<ValueType> &a_half, lama::DenseVector<ValueType> &b_half, lama::DenseVector<ValueType> &kInv_half,
                                                                                      IndexType i)
{

    utilskernel::LArray<ValueType> *k_LA = &kInv.getLocalValues();
    hmemo::WriteAccess<ValueType> write_k(*k_LA);

    utilskernel::LArray<ValueType> *b_LA = &b.getLocalValues();
    hmemo::WriteAccess<ValueType> write_b(*b_LA);

    utilskernel::LArray<ValueType> *a_LA = &a.getLocalValues();
    hmemo::WriteAccess<ValueType> write_a(*a_LA);

    utilskernel::LArray<ValueType> *k_half_LA = &kInv_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_k_half(*k_half_LA);

    utilskernel::LArray<ValueType> *b_half_LA = &b_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_b_half(*b_half_LA);

    utilskernel::LArray<ValueType> *a_half_LA = &a_half.getLocalValues();
    hmemo::WriteAccess<ValueType> write_a_half(*a_half_LA);

    write_a[i] = 0.0;
    write_b[i] = 0.0;
    write_k[i] = 1.0;

    write_a_half[i] = 0.0;
    write_b_half[i] = 0.0;
    write_k_half[i] = 1.0;

    write_k.release();
    write_a.release();
    write_b.release();
    write_k_half.release();
    write_a_half.release();
    write_b_half.release();
}

/*! \brief Application of the CPML
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param Vec DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 \param a CPML coefficient a
 \param b CPML coefficient b
 \param kInv reciprocal of CPML coefficient k
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::applyCPML(lama::Vector &Vec, lama::DenseVector<ValueType> &Psi, lama::DenseVector<ValueType> &a, lama::DenseVector<ValueType> &b, lama::DenseVector<ValueType> &kInv)
{

    update_PmlTemp = Vec;
    Psi.scale(b);
    Psi += update_PmlTemp.scale(a);
    Vec.scale(kInv);
    Vec += Psi;
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML<float>;
