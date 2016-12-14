
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "../../Acquisition/Acquisition.hpp"
#include "../../Acquisition/Seismogram.hpp"
#include "../../Acquisition/SeismogramHandler.hpp"

#include "../../Modelparameter/Modelparameter.hpp"
#include "../../Wavefields/Wavefields.hpp"

#include "FDTDacoustic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        namespace SourceReceiverImpl {
            
            template<typename ValueType>
            class FDTD2Dacoustic : public FDTDacoustic<ValueType>
            {
            public:
                
                FDTD2Dacoustic()=delete;
                ~FDTD2Dacoustic(){};
                
                using FDTDacoustic<ValueType>::FDTDacoustic;

            };
            
        }
    }
}


