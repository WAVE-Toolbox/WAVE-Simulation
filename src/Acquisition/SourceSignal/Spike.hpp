#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "SourceSignal.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        namespace SourceSignal {
            
            /*! \brief Class to create a Spike sourcesignals
             *
             */
            template <typename ValueType>
            class Spike : public SourceSignal<ValueType> {
            public:
                
                explicit Spike(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift);
                
                ///! Destructor
                ~Spike(){};
                
                void calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) override;
                
            };
        }
    }
}

/*! \brief Constructor generating a Spike signal
 *
 \param signal Allocated vector to store Spike signal
 \param NT Number of time steps
 \param DT Temporal time step interval
 \param FC Will not be used
 \param AMP Amplitude
 \param Tshift Time to shift wavelet
 */
template<typename ValueType>
KITGPI::Acquisition::SourceSignal::Spike<ValueType>::Spike(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift ){
    calc(signal, NT, DT, FC, AMP, Tshift );
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
template<typename ValueType>
void KITGPI::Acquisition::SourceSignal::Spike<ValueType>::calc(lama::DenseVector<ValueType>& signal, IndexType NT, ValueType DT, ValueType /*FC*/, ValueType AMP, ValueType Tshift )
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
