#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <cmath>
#include <valarray>

namespace KITGPI {
    
    namespace Acquisition {
         
        /*! \brief Abstract class to create sourcesignals
         *
         * This is class holds several methods to generate different sourcesignals.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class Sourcesignal
        {
            
        protected:
            void Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
	    void FGaussian(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
	    void Spike(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType AMP, ValueType Tshift);
	    void sinthree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
	    void intgsinthree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
	    void sinw(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
        };
    }
}


/*! \brief Generating a Ricker signal
 *
 \param signal Allocated vector to store Ricker signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
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

/*! \brief Generating a First derivative of a Gaussian (FGaussian)
 *
 \param signal Allocated vector to store FGaussian signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::FGaussian(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
{
    
    /*
     *  t=0:DT:(NT*DT-DT);
     *  tau=pi*FC*(t-1.2/FC);
     *  signal=AMP*2*tau.*exp(-tau*tau);
     */
    lama::DenseVector<ValueType> t(NT, ValueType(0), DT);
    lama::DenseVector<ValueType> help( t.size(), 1.2 / FC +Tshift);
    lama::DenseVector<ValueType> tau( t - help );
    tau *= M_PI * FC ;
    
    /* this is for source[i] = AMP * (-2) * tau[i].*exp(-tau[i] * tau[i]); */

    help = -2.0 * tau;
    tau = -1.0* tau * tau;
    tau.exp();
    signal = lama::Scalar(AMP) * help * tau;
    
}


/*! \brief Generating a Spike/Delta wavelet
 *
 \param signal Allocated vector to store a Spike signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::Spike(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType AMP, ValueType Tshift )
{
    
    /*
     *  Spike;
     */
    scai::lama::Scalar temp_spike;
    IndexType time_index;
    lama::DenseVector<ValueType> help( NT, 0.0);
    
    /* this is for source[i] = 1.0 when t=tshift/dt; */
    temp_spike=1.0;
    time_index=floor(Tshift/DT);
    help.setValue(time_index,temp_spike);

    signal = lama::Scalar(AMP) * help;
    
}

/*! \brief Generating a sinus raised to the power of three (sinthree)
 *
 \param signal Allocated vector to store sinthree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::sinthree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
{
    
     /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=pi*FC*(t-Tshift);
     *  signal=(sin(tau))^3;
     */

    lama::DenseVector<ValueType> zero( NT, 0.0 );
    
    double temp;
    IndexType time_index1,time_index2,i,count;
    
    time_index1 = floor(Tshift/DT);
    time_index2 = time_index1 + floor(1.0/FC/DT);
    
    
    /* this is for source[i] = (sin(PI*(t-Tshift)*FC))^3 when t>=tshift && t<=tshift+1.0/FC; */
    count=0;
    for (i=time_index1; i<=time_index2; i++) {
	  temp=count * DT * M_PI * FC ;
	  temp=sin(temp);
	  temp=pow(temp,3);
          zero.setValue(i,temp);
	  count++;
    }

    signal = lama::Scalar(AMP) * zero;
    
}

/*! \brief Generating a wavelet composed by integral of sinus raised to the power of three (sinthree)
 *
 \param signal Allocated vector to store integral sinthree signal (intgsinthree)
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::intgsinthree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
{
    
     /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=pi*FC*(t-Tshift);
     *  signal=(0.5-0.75*cos(tau)+0.25*(sin(tau))^3)/FC/0.75/pi;
     */

    lama::DenseVector<ValueType> zero( NT, 0.0 );
    lama::DenseVector<ValueType> help( NT, 0.0 );
    lama::DenseVector<ValueType> half( NT, 0.5 );
    
    double temp;
    IndexType time_index1,time_index2,i,count;
    
    time_index1 = floor(Tshift/DT);
    time_index2 = time_index1 + floor(1.0/FC/DT);
    
    
    /* this is for source[i] = (0.5-0.75*cos(tau[i])+0.25*(cos(tau[i]))^3)/FC/0.75/pi when t>=tshift && t<=tshift+1.0/FC; */
    count=0;
    for (i=time_index1; i<=time_index2; i++) {
	  temp=count * DT * M_PI * FC ;
	  temp=cos(temp);
	  help.setValue(i,temp);
	  temp=pow(temp,3);
          zero.setValue(i,temp);
	  count++;	  
    }
    help = 0.75 *help;
    zero = 0.25 * zero;
    zero = zero - help;
    zero = half + zero;
    zero = zero/FC;
    zero = zero/M_PI;
    zero = zero/0.75;

    signal = lama::Scalar(AMP) * zero;
    
}

/*! \brief Generating a wavelet composed by a combination of two sin functions (sinw)
 *
 \param signal Allocated vector to to store wavelet signal (sinw)
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::Sourcesignal<ValueType>::sinw(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
{
    
     /*
     *  t=0:DT:(NT*DT-DT);
     *  when t>=tshift && t<=tshift+1.0/FC;
     *  tau=2*pi*FC*(t-Tshift);
     *  signal=sin(tau)-0.5*sin(2*tau);
     */

    lama::DenseVector<ValueType> zero( NT, 0.0 );
    lama::DenseVector<ValueType> help( NT, 0.0 );
    lama::DenseVector<ValueType> half( NT, 0.5 );
    
    double temp,temp2;
    IndexType time_index1,time_index2,i,count;
    
    time_index1 = floor(Tshift/DT);
    time_index2 = time_index1 + floor(1.0/FC/DT);
    
    
    /* this is for source[i] = sin(tau[i])-0.5*sin(2*tau[i]) when t>=tshift && t<=tshift+1.0/FC; */
    count=0;
    for (i=time_index1; i<=time_index2; i++) {
	  temp= 2.0 * count * DT * M_PI * FC ;
	  temp2=sin(temp);
	  help.setValue(i,temp2);
	  temp = 2.0 * temp;
	  temp2 = sin(temp);
          zero.setValue(i,temp2);
	  count++;	  
    }
    zero = zero/2.0;
    help = help - zero;

    signal = lama::Scalar(AMP) * help;
    
}

