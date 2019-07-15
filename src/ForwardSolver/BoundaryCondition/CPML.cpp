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
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::initVector(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr const ctx, scai::dmemo::DistributionPtr const dist)
{
    vector.setContextPtr(ctx);
    vector.setSameValue(dist, 0.0);
}

/*! \brief set CPML coefficients
 * 
 * method to set cpml coefficients for a given gridpoint
 *  
 \f{eqnarray*}{
        a =& \frac{d}{d+\alpha} (b-1) \\
        b =& e^{-\Delta t (d+\alpha)} \\
        d =& d_0 (\frac{w_b-g_d}{w_b})^N \\
        d_0 =& \frac{-(N+1)v_{max} \ln (8\cdot 10^{-4})}{2 w_b \Delta x} \\
        \alpha =& \pi f_c (1-(\frac{w_b-g_d}{w_b})^N)
 \f}
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::calcCoeffCPML(std::vector<ValueType> &a, std::vector<ValueType> &b, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, ValueType const DT, ValueType const DH, bool const shiftGrid)
{
    IndexType BoundaryWidth = a.size();

    ValueType shift = 0;
    if (shiftGrid) {
        shift = 0.5;
    }

    ValueType RCoef = 0.0008;

    ValueType alpha_max_Pml = 2.0 * M_PI * (CenterFrequencyCPML / 2.0);

    ValueType d0 = -(NPower + 1) * VMaxCPML * log(RCoef) / (2.0 * BoundaryWidth * DH);

    ValueType PositionNorm, alpha_prime, d = 0.0;

    for (IndexType i = 0; i < BoundaryWidth; i++) {
        PositionNorm = (ValueType)(BoundaryWidth - i - shift) / BoundaryWidth;
        d = d0 * pow(PositionNorm, NPower);
        alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
        b[i] = exp(-(d + alpha_prime) * DT);
        /* avoid division by zero outside the PML */
        if (std::abs(d) > 1.0e-6) {
            a[i] = d * (b[i] - 1.0) / (d + alpha_prime);
        } else {
            a[i] = 0.0;
        }
    }
}

/*! \brief Application of the CPML. Replace FD-operators \f$\partial\f$:
  \f{eqnarray*}{
        \bar{\partial} = \partial + \Psi \\
        \Psi^n = b\Psi^{n-1} + a \partial^{n+\frac{1}{2}} 
 \f} 
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param Vec DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 \param a CPML coefficient a
 \param b CPML coefficient b
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::applyCPML(lama::DenseVector<ValueType> &Vec, VectorType &Psi, VectorType const &a, VectorType const &b)
{
    temp = a;
    Psi *= b;
    temp *= Vec;
    Psi += temp;
    Vec += Psi;
}

template class KITGPI::ForwardSolver::BoundaryCondition::CPML<double>;
template class KITGPI::ForwardSolver::BoundaryCondition::CPML<float>;
