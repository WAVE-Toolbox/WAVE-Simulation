#include "CPML.hpp"
using namespace scai;

/*! \brief Reset a single Vector to zero.
*/
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::resetVector(scai::lama::Vector<ValueType> &vector)
{
    vector = 0.0;
}

/*! \brief Intitialisation of a single vector.
*
* This method will set the context, allocate the the wavefield and set the field to zero.
*/
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::initVector(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);

    resetVector(vector);
}

/*! \brief set CPML coefficients
 * 
 * method to set cpml coefficients for a given gridpoint
 *  
 \f{eqnarray*}{
        a =& \frac{d}{k(d+k\alpha)} (b-1) \\
        b =& e^{-\Delta t (d/k+ \alpha)} \\
        d =& d_0 (\frac{w_b-g_d}{w_b})^N \\
        d_0 =& \frac{-(N+1)v_{max} \ln (8\cdot 10^{-4})}{2 w_b \Delta x} \\
        k =& 1 + (k_{max} -1) (\frac{w_b-g_d}{w_b})^N \\
        \alpha =& \pi f_c (1-(\frac{w_b-g_d}{w_b})^N)
 \f}
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

    auto write_k = hmemo::hostWriteAccess(kInv.getLocalValues());
    auto write_b = hmemo::hostWriteAccess(b.getLocalValues());
    auto write_a = hmemo::hostWriteAccess(a.getLocalValues());
    auto write_k_half = hmemo::hostWriteAccess(kInv_half.getLocalValues());
    auto write_b_half = hmemo::hostWriteAccess(b_half.getLocalValues());
    auto write_a_half = hmemo::hostWriteAccess(a_half.getLocalValues());

    /* left boundary */
    if (coord < BoundaryWidth) {
        PositionNorm = (ValueType)(BoundaryWidth - gdist) / BoundaryWidth;
        d = d0 * pow(PositionNorm, NPower);
        k_temp = 1.0 + (KMaxCPML - 1.0) * pow(PositionNorm, NPower);
        alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
        b_temp = exp(-(d / k_temp + alpha_prime) * DT);
        /* avoid division by zero outside the PML */
        if (std::abs(d) > 1.0e-6) {
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
        if (std::abs(d) > 1.0e-6) {
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
    if (std::abs(d) > 1.0e-6) {
        a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));
    } else
        a_temp = 0.0;
    write_a_half[i] = a_temp;
    write_b_half[i] = b_temp;
    write_k_half[i] = 1 / k_temp;
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

    auto write_k = hmemo::hostWriteAccess(kInv.getLocalValues());
    auto write_b = hmemo::hostWriteAccess(b.getLocalValues());
    auto write_a = hmemo::hostWriteAccess(a.getLocalValues());
    auto write_k_half = hmemo::hostWriteAccess(kInv_half.getLocalValues());
    auto write_b_half = hmemo::hostWriteAccess(b_half.getLocalValues());
    auto write_a_half = hmemo::hostWriteAccess(a_half.getLocalValues());

    write_a[i] = 0.0;
    write_b[i] = 0.0;
    write_k[i] = 1.0;

    write_a_half[i] = 0.0;
    write_b_half[i] = 0.0;
    write_k_half[i] = 1.0;
}

/*! \brief Application of the CPML. Replace FD-operators \f$\partial\f$:
  \f{eqnarray*}{
        \bar{\partial} = \frac{1}{k} \partial + \Psi \\
        \Psi^n = b\Psi^{n-1} + a \partial^{n+\frac{1}{2}} 
 \f} 
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
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::applyCPML(lama::Vector<ValueType> &Vec, lama::Vector<ValueType> &Psi, lama::Vector<ValueType> &a, lama::Vector<ValueType> &b, lama::Vector<ValueType> &kInv)
{
    // the following directive guarantees that two sparse vectors with same number of entries have the same pattern

    SCAI_SPARSE_VECTOR_SAME_PATTERN

    temp = a;
    Psi *= b;
    temp *= Vec;
    Psi += temp;
    Vec *= kInv;
    Vec += Psi;
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML<float>;
