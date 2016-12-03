#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        namespace SourceSignal {
            
            /*! \brief Class to create a Ricker sourcesignals
             *
             */
            template <typename ValueType>
            class Ricker : public SourceSignal<ValueType> {
            public:
                
                explicit Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
                
                ///! Destructor
                ~Ricker(){};
                
                void calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
                
            };
        }
    }
}

/*! \brief Constructor generating a Ricker signal
 *
 \param signal Allocated vector to store Ricker signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
KITGPI::Acquisition::SourceSignal::Ricker<ValueType>::Ricker(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift ){
    calc(signal, NT, DT, FC, AMP, Tshift );
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
void KITGPI::Acquisition::SourceSignal::Ricker<ValueType>::calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
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
