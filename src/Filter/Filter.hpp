#pragma once
#include <scai/common/Complex.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/fft.hpp>

#include <algorithm>
#include <string>

#include "../Common/Common.hpp"

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

            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;

            void init(ValueType dt, scai::IndexType nt);
            void calc(std::string transFcnFmly, std::string filterType, scai::IndexType order, ValueType fc1, ValueType fc2 = 0.0);
            void apply(scai::lama::DenseVector<ValueType> &signal) const;
            void apply(scai::lama::DenseMatrix<ValueType> &signal) const;
            //
          private:
            scai::IndexType zeroPadding;
            scai::lama::DenseVector<ComplexValueType> transFcn;
            ValueType df;
            ValueType fNyquist;

            scai::IndexType calcZeroPadding(scai::IndexType nt);
            void calcFrequencyVector(scai::lama::DenseVector<ValueType> &frequencyVector);
            void calcFrequencyMat(std::string filterType, scai::IndexType order, ValueType fc, scai::lama::DenseVector<ValueType> &frequencyVector, scai::lama::DenseMatrix<ComplexValueType> &frequencyMat);

            void calcButterPoly(scai::IndexType order, scai::lama::DenseVector<ValueType> &poly);

            void calcButterworthFilt(std::string filterType, scai::IndexType order, ValueType fc1, ValueType fc2 = 0.0);
            void calcButterworthHp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc);
            void calcButterworthLp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc);
            void calcButterworthBp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc1, ValueType fc2);
        };
    }
}
