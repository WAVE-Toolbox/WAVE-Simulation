#include "Sin.hpp"
using namespace scai;

/*! \brief Constructor generating a Sin signal
 *
 \param signal Allocated vector to store Sin signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
KITGPI::Acquisition::SourceSignal::Sin<ValueType>::Sin(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{
    calc(signal, NT, DT, FC, AMP, Tshift);
}

/*! \brief Generating a Sin signal
 *
 \param signal Allocated vector to store Sin signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template <typename ValueType>
void KITGPI::Acquisition::SourceSignal::Sin<ValueType>::calc(scai::lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(DT > 0, "DT is < 0: No valid argument!");
    SCAI_ASSERT_ERROR(FC > 0, "FC is < 0: No valid argument!");

    /*
     *  t=0:DT:(NT*DT-DT);
     *  tau=2*pi*FC*(t-Tshift);
     *  signal=sin(tau);
     */

    lama::DenseVector<ValueType> tau = lama::linearDenseVector<ValueType>(NT, 0, DT);
    tau -= Tshift;
    tau *= 2.0 * M_PI * FC;
    tau = scai::lama::sin(tau);
    
    signal = AMP * tau;
}

template class KITGPI::Acquisition::SourceSignal::Sin<float>;
template class KITGPI::Acquisition::SourceSignal::Sin<double>;
