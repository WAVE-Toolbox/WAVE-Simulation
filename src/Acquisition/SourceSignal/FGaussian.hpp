#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        namespace SourceSignal {
            
            /*! \brief Class to create a FGaussian sourcesignals
             *
             */
            template <typename ValueType>
            class FGaussian : public SourceSignal<ValueType> {
            public:
                
                FGaussian(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
                
                ///! Destructor
                ~FGaussian(){};
                
                void calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
                
            };
        }
    }
}

/*! \brief Constructor generating a FGaussian signal
 *
 \param signal Allocated vector to store IntgSinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
KITGPI::Acquisition::SourceSignal::FGaussian<ValueType>::FGaussian(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift ){
    calc(signal, NT, DT, FC, AMP, Tshift );
}

/*! \brief Generating a FGaussian signal
 *
 \param signal Allocated vector to store IntgSinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
void KITGPI::Acquisition::SourceSignal::FGaussian<ValueType>::calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
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
