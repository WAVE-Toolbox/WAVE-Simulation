
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#pragma once

/*! \brief Abstract class to create sourcesignals
 * This is class holds several methods to generate different sourcesignals.
 * As this class is an abstract class, all methods are protected.
 */
template <typename ValueType>
class Sourcesignal
{

protected:
    
    void Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
};


/*! \brief Generating a Ricker signal
 *
 \param signal Allocated vector to to store Ricker signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void Sourcesignal<ValueType>::Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
{
    
    /*
     *  t=0:DT:(NT*DT-DT);
     *  tau=pi*FC*(t-1.5/FC);
     *  signal=AMP*(1-2*tau.^2).*exp(-tau.^2);
     */
    lama::DenseVector<ValueType> t(NT, ValueType(0), DT);
    lama::DenseVector<ValueType> help( t.size(), 1.5 / FC +Tshift);
    lama::DenseVector<ValueType> tau( t - help );
    tau *= M_PI * FC;
    
    /* this is for source[i] = AMP * ( 1.0 - 2.0 * tau[i] * tau[i] * exp( -tau[i] * tau[i] ) ); */
    lama::DenseVector<ValueType> one( signal.size(), 1.0 );
    help = tau * tau;
    tau = -1.0 * help;
    tau.exp();
    help = one - 2.0 * help;
    signal = lama::Scalar(AMP) * help * tau;
    
}
