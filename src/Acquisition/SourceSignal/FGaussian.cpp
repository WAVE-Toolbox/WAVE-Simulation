#include "FGaussian.hpp"
using namespace scai;

/*! \brief Constructor generating a FGaussian signal
 *
 \param signal Allocated vector to store IntgSinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::FGaussian<ValueType>::FGaussian(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a FGaussian signal
 *
 \param signal Allocated vector to store FGaussian signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::FGaussian<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "FC is < 0: No valid argument!");

    /*
     *  t=0:DT:(NT*DT-DT);
     *  tau=pi*FC*(t-1.2/FC-Tshift);
     *  signal=-AMP*2*tau.*exp(-tau*tau);
     */
    auto t = lama::linearDenseVector<ValueType>(NT, 0, DT);
    auto help = lama::fill<lama::DenseVector<ValueType>>(t.size(), 1.2 / FC + Tshift);
    auto tau = lama::eval<lama::DenseVector<ValueType>>(t - help);
    tau *= M_PI * FC;

    /* this is for source[i] = AMP * (-2) * tau[i].*exp(-tau[i] * tau[i]); */

    help = -2.0 * tau;
    tau = -1.0 * tau * tau;
    tau = exp(tau);
    signal = AMP * help * tau;
}

template class KITGPI::Acquisition::SourceSignal::FGaussian<float>;
template class KITGPI::Acquisition::SourceSignal::FGaussian<double>;
