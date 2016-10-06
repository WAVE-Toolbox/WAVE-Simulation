
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "../Derivatives/Derivatives.hpp"

template<typename ValueType>
class ForwardSolver
{
public:
    
    ForwardSolver(){};
    ~ForwardSolver(){};
    
    virtual void run(Receivers<ValueType>& receiver, Sources<ValueType>& sources, Modelparameter<ValueType>& model, Wavefields<ValueType>& wavefield, Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm)=0;
            
};
