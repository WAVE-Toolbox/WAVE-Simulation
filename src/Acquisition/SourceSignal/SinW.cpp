#include "SinW.hpp"
using namespace scai;

/*! \brief Constructor generating a SinW signal
 *
 \param signal Allocated vector to store SinW signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::SinW<ValueType>::SinW(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a SinW signal
 *
 \param signal Allocated vector to store SinW signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::SinW<ValueType>::calc(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "DT is < 0: No valid argument!");

    /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=2*pi*FC*(t-Tshift);
     *  signal=sin(tau)-0.5*sin(2*tau);
     */

    lama::DenseVector<ValueType> zero(NT, 0.0);
    lama::DenseVector<ValueType> help(NT, 0.0);
    lama::DenseVector<ValueType> half(NT, 0.5);

    double temp, temp2;
    IndexType time_index1, time_index2, i, count;

    time_index1 = floor(Tshift / DT);
    time_index2 = time_index1 + floor(1.0 / FC / DT);

    /* this is for source[i] = sin(tau[i])-0.5*sin(2*tau[i]) when t>=tshift && t<=tshift+1.0/FC; */
    count = 0;
    for (i = time_index1; i <= time_index2; i++) {
        temp = 2.0 * count * DT * M_PI * FC;
        temp2 = sin(temp);
        help.setValue(i, temp2);
        temp = 2.0 * temp;
        temp2 = sin(temp);
        zero.setValue(i, temp2);
        count++;
    }
    zero = zero / 2.0;
    help = help - zero;

    signal = lama::Scalar(AMP) * help;
}

template class KITGPI::Acquisition::SourceSignal::SinW<float>;
template class KITGPI::Acquisition::SourceSignal::SinW<double>;
