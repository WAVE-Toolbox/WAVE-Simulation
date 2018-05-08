#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/lama/fft.hpp>
#include <scai/common/Complex.hpp>

#include <algorithm>
#include <string> 

namespace KITGPI
{

    //! \brief Filter namespace
    namespace Filter
    {
        
        //! \brief Class to handle frequency filtering.
        template <typename ValueType>
        class Filter
        {
        public:
            //! Default constructor
            Filter(){};
            //! Default destructor
            ~Filter(){};
            
            void init(ValueType dt, scai::IndexType nt);
            void calc(std::string transFcnFmly, std::string filterType, ValueType fc, scai::IndexType order);
            void apply(scai::lama::DenseVector<ValueType> &signal);
        private:
            scai::IndexType zeroPadding;
            scai::lama::DenseVector<ValueType> transFcn;
            ValueType df;
            ValueType fNyquist;
            
            scai::IndexType calcZeroPadding(scai::IndexType nt);
            void calcFrequencyVector(scai::lama::DenseVector<ValueType> &frequencyVector);
            
            void calcButterworthFilt(std::string filterType, ValueType fc, scai::IndexType order);
            
        };
    }
}