#include "Ricker_GprMax.hpp"
using namespace scai;

/*! \brief Constructor generating a Ricker_GprMax signal
 *
 \param signal Allocated vector to store Ricker_GprMax signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::Ricker_GprMax<ValueType>::Ricker_GprMax(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a Ricker_GprMax signal
 *
 \param signal Allocated vector to store Ricker_GprMax signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::Ricker_GprMax<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "DT is < 0: No valid argument!");

    /*
     *  t=0:DT:(NT*DT-DT);
     *  zeta=2*PI*PI*FC*FC;
     *  tau=t-(1/FC+Tshift);
     *  signal=AMP*(-2)*zeta*(exp(1/(2*zeta))).^(1/2)*exp(-zeta*tau.^2)*tau;
     */
    auto t = lama::linearDenseVector<ValueType>(NT, 0, DT);
    lama::DenseVector<ValueType> help(t.size(), 1 / FC + Tshift);
    auto tau = lama::eval<lama::DenseVector<ValueType>>(t - help);
    lama::DenseVector<ValueType> zeta(t.size(), 2.0 * M_PI * M_PI * FC * FC);


    help = tau * tau;
    help = zeta * help;
    help = -1.0 * help;
    help = exp(help);
    tau = help * tau;
    help = 1.0 / zeta;
    help = help / 4;
    help = exp(help);
    zeta = zeta *help;
    zeta = -2.0 * zeta;
    signal = AMP * zeta * tau;
}

template class KITGPI::Acquisition::SourceSignal::Ricker_GprMax<float>;
template class KITGPI::Acquisition::SourceSignal::Ricker_GprMax<double>;
