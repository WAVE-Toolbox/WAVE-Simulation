#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

namespace KITGPI
{

    namespace Acquisition
    {

        /*! SourceSignal class for calculation of synthetic wavelets */
        namespace SourceSignal
        {

            /*! \brief Abstract class to create SourceSignals
             *
             * As this class is an abstract class.
             */
            template <typename ValueType>
            class SourceSignal
            {
              public:
                //! Constructor
                SourceSignal(){};

                //! Destructor
                ~SourceSignal(){};

                /*! \brief Generating a source signal
                 *
                 \param signal Allocated vector to store source signal
                 \param NT Number of time steps
                 \param DT Temporal time step interval
                 \param FC Central frequency
                 \param AMP Amplitude
                 \param Tshift Time to shift wavelet
                 */
                virtual void calc(lama::DenseVector<ValueType> &signal, IndexType NT, ValueType DT, ValueType FC, ValueType AMP, ValueType Tshift) = 0;
            };
        }
    }
}
