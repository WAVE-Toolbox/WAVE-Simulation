#include "SinThree.hpp"
using namespace scai;

/*! \brief Constructor generating a SinThree signal
 *
 \param signal Allocated vector to store SinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::SinThree<ValueType>::SinThree(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a SinThree signal
 *
 \param signal Allocated vector to store SinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::SinThree<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "FC is < 0: No valid argument!");

    /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=pi*FC*(t-Tshift);
     *  signal=(sin(tau))^3;
     */

    lama::DenseVector<ValueType> zero(NT, 0.0);

    double temp;
    IndexType time_index1, time_index2, i, count;

    time_index1 = floor(Tshift / DT);
    time_index2 = time_index1 + floor(1.0 / FC / DT);
    
    SCAI_ASSERT(time_index1<signal.size(),"Signal is shifted to Tshift > Tsimulation");
    if (time_index2>signal.size()){
        time_index2=signal.size()-1;
    }

    /* this is for source[i] = (sin(PI*(t-Tshift)*FC))^3 when t>=tshift && t<=tshift+1.0/FC; */
    count = 0;
    for (i = time_index1; i <= time_index2; i++) {
        temp = count * DT * M_PI * FC;
        temp = sin(temp);
        temp = pow(temp, 3);
        zero.setValue(i, temp);
        count++;
    }

    signal = AMP * zero;
}

template class KITGPI::Acquisition::SourceSignal::SinThree<float>;
template class KITGPI::Acquisition::SourceSignal::SinThree<double>;
