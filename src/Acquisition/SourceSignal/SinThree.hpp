#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        namespace SourceSignal {
            
            /*! \brief Class to create a SinThree sourcesignals
             *
             */
            template <typename ValueType>
            class SinThree : public SourceSignal<ValueType> {
            public:
                
                explicit SinThree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
                
                ///! Destructor
                ~SinThree(){};
                
                void calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
                
            };
        }
    }
}

/*! \brief Constructor generating a SinThree signal
 *
 \param signal Allocated vector to store SinThree signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Central frequency
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
KITGPI::Acquisition::SourceSignal::SinThree<ValueType>::SinThree(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift ){
    calc(signal, NT, DT, FC, AMP, Tshift );
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
template<typename ValueType>
void KITGPI::Acquisition::SourceSignal::SinThree<ValueType>::calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift )
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
        temp=count * DT * M_PI *  FC ;
        temp=sin(temp);
        temp=pow(temp,3);
        zero.setValue(i,temp);
        count++;
    }
    
    signal = lama::Scalar(AMP) * zero;
    
    
}
