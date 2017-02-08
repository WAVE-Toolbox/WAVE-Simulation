#include "IntgSinThree.hpp"

/*! \brief Constructor generating a IntgSinThree signal
 *
 \param signal Allocated vector to store IntgSinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::IntgSinThree<ValueType>::IntgSinThree(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a IntgSinThree signal
 *
 \param signal Allocated vector to store IntgSinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::IntgSinThree<ValueType>::calc(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "DT is < 0: No valid argument!");
    
    /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=pi*FC*(t-Tshift);
     *  signal=(0.5-0.75*cos(tau)+0.25*(sin(tau))^3)/FC/0.75/pi;
     */

    lama::DenseVector<ValueType> zero(NT, 0.0);
    lama::DenseVector<ValueType> help(NT, 0.0);
    lama::DenseVector<ValueType> half(NT, 0.5);

    double temp;
    IndexType time_index1, time_index2, i, count;

    time_index1 = floor(Tshift / DT);
    time_index2 = time_index1 + floor(1.0 / FC / DT);

    /* this is for source[i] = (0.5-0.75*cos(tau[i])+0.25*(cos(tau[i]))^3)/FC/0.75/pi when t>=tshift && t<=tshift+1.0/FC; */
    count = 0;
    for (i = time_index1; i <= time_index2; i++) {
        temp = count * DT * M_PI * FC;
        temp = cos(temp);
        help.setValue(i, temp);
        temp = pow(temp, 3);
        zero.setValue(i, temp);
        count++;
    }
    help = 0.75 * help;
    zero = 0.25 * zero;
    zero = zero - help;
    zero = half + zero;
    zero = zero / FC;
    zero = zero / M_PI;
    zero = zero / 0.75;

    signal = lama::Scalar(AMP) * zero;
}

template class KITGPI::Acquisition::SourceSignal::IntgSinThree<float>;
template class KITGPI::Acquisition::SourceSignal::IntgSinThree<double>;
