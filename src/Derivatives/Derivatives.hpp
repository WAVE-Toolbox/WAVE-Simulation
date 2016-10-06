
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

namespace KITGPI {
    
    //! \brief Derivatives namespace
    namespace Derivatives {
        
        template<typename ValueType>
        class Derivatives
        {
        public:
            
            Derivatives(){};
            ~Derivatives(){};
            
            virtual lama::CSRSparseMatrix<ValueType>* getA()=0;
            virtual lama::CSRSparseMatrix<ValueType>* getB()=0;
            virtual lama::CSRSparseMatrix<ValueType>* getC()=0;
            virtual lama::CSRSparseMatrix<ValueType>* getD()=0;
            virtual lama::CSRSparseMatrix<ValueType>* getE()=0;
            virtual lama::CSRSparseMatrix<ValueType>* getF()=0;
            
            virtual void initializeMatrices(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, dmemo::CommunicatorPtr comm )=0;
            
            
        };
    }
}
