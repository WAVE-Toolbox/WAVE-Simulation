#include "Spike.hpp"
using namespace scai;

/*! \brief Constructor generating a Spike signal
 *
 \param signal Allocated vector to store Spike signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Will not be used
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::Spike<ValueType>::Spike(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a Spike signal
 *
 \param signal Allocated vector to store Spike signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Will not be used
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::Spike<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType /*FC*/, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");

    /*
     *  Spike;
     */
    scai::lama::Scalar temp_spike;
    IndexType time_index;
    lama::DenseVector<ValueType> help(NT, 0.0);

    /* this is for source[i] = 1.0 when t=tshift/dt; */
    temp_spike = 1.0;
    time_index = floor(Tshift / DT);
    help.setValue(time_index, temp_spike);

    signal = lama::Scalar(AMP) * help;
}

template class KITGPI::Acquisition::SourceSignal::Spike<double>;
template class KITGPI::Acquisition::SourceSignal::Spike<float>;
