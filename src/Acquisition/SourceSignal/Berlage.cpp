#include "Berlage.hpp"
using namespace scai;

/*! \brief Constructor generating a Berlage signal
 *
 \param signal Allocated vector to store Berlage signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::Berlage<ValueType>::Berlage(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a Berlage signal
 *
 \param signal Allocated vector to store Berlage signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::Berlage<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "FC is < 0: No valid argument!");
    
    /*
     *  Aldridge, David F. "The Berlage wavelet." Geophysics 55.11 (1990): 1508-1511.
     *  t=0:DT:(NT*DT-DT); t=t-1.0/FC-Tshift;
     *  tau=2*pi*FC*t; n=2; alpha=2*FC;     
     *  signal=cos(tau);
     *  tau=heaviside(t).*t.^n.*exp(-alpha.*t);
     *  signal*=tau;
     *  signal*=AMP/signal.maxNorm();
     */
    auto t = lama::linearDenseVector<ValueType>(NT, 0, DT);
    lama::DenseVector<ValueType> help(t.size(), 1.0 / FC + Tshift);
    auto tau = lama::eval<lama::DenseVector<ValueType>>(t - help);
    ValueType n = 2;
    ValueType alpha = 2 * FC;
    
    t = 2 * M_PI * FC * tau;
    signal = lama::cos(t);

    lama::DenseVector<ValueType> heaviside;
    heaviside.unaryOp(tau, common::UnaryOp::SIGN);
    heaviside += 1;
    heaviside.unaryOp(heaviside, common::UnaryOp::SIGN);
    help = lama::pow(tau, n);
    tau *= -alpha;
    tau = lama::exp(tau);
    signal *= help;
    signal *= tau;
    signal *= heaviside;
    signal *= AMP / signal.maxNorm();
}

template class KITGPI::Acquisition::SourceSignal::Berlage<float>;
template class KITGPI::Acquisition::SourceSignal::Berlage<double>;
