#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
// #include <scai/common/Stencil.hpp>
// #include <scai/lama.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
// #include <scai/lama/matrix/StencilMatrix.hpp>
// #include <scai/lama/matrix/HybridMatrix.hpp>
#include <scai/tracing.hpp>

#include <cmath>

namespace KITGPI
{
    //! \brief Common namespace
    namespace Common
    {

        /*! \brief Searches for all values in searchVector which are related to threshold by relation compareType and replaces them with replaceValue.
        *
        \param searchVector The vector to be searched. 
        \param threshold The threshhold which the values are compared to.
        \param replaceValue The value by which the enties found in searchVector are replaced with.
        \param compareType The relation the values in searchVector and threshold should have. Possible values are 1 := <, 2 := >, 3 := <=, 4 := >=, 5 := ==
        */
        template <typename ValueType>
        void searchAndReplace(scai::lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, scai::IndexType compareType)
        {

            // needs rework with new LAMA features
            scai::hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);

            switch (compareType) {
            case 1: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] < threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 2: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] > threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 3: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] <= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 4: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] >= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 5: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] == threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            default:
                COMMON_THROWEXCEPTION("Invalid compareType. Has to be < 6 but is " << compareType)
            }

            write_searchVector.release();
        }

        /*! \brief Replaces NaN and Inf values in a vector by a given value
        *
        \param searchVector Input vector
        \param replaceValue Value NaN and Inf are set to
        */
        template <typename ValueType>
        void replaceInvalid(scai::lama::DenseVector<ValueType> &searchVector, ValueType replaceValue)
        {

            // needs rework with new LAMA features
            scai::hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);

            for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (std::isnan(write_searchVector[i]) || write_searchVector[i] == std::numeric_limits<ValueType>::infinity() || -write_searchVector[i] == std::numeric_limits<ValueType>::infinity()) { //std::isinf doesn't work for whatever reason
                    write_searchVector[i] = replaceValue;
                }
            }

            write_searchVector.release();
        }

        /*! \brief Calculate the next power of two.
        \param nt number the next power of two should be calculated for
        */
        template <typename ValueType>
        scai::IndexType calcNextPowTwo(scai::IndexType nt)
        {
            ValueType temp = scai::common::Math::log(ValueType(nt));
            temp /= scai::common::Math::log(2.0);
            temp = scai::common::Math::ceil(temp);
            temp = scai::common::Math::pow(ValueType(2.0), temp);
            return temp;
        }

        /*! \brief Calculate a matrix which resamples the columns.
        \param rMat resampling matrix
        \param numCols number of samples in one row
        \param resamplingCoeff resampling coefficient
        */
        template <typename ValueType>
        void calcResampleMat(scai::lama::CSRSparseMatrix<ValueType> &rMat, scai::IndexType numCols, ValueType resamplingCoeff)
        {

            scai::lama::MatrixAssembly<ValueType> assembly;

            scai::IndexType numColsNew = scai::IndexType(scai::common::Math::floor<ValueType>(ValueType(numCols - 1) / ValueType(resamplingCoeff))) + 1; // number of samples after resampling

            scai::IndexType columnIndex = 0;
            ValueType value = 0.0;
            ValueType sampleCoeff = 1.0;
            for (scai::IndexType rowIndex = 0; rowIndex < numColsNew; rowIndex++) {

                ValueType relativeIndex = rowIndex * resamplingCoeff;
                //leftValue
                columnIndex = scai::common::Math::floor<ValueType>(relativeIndex);
                if (columnIndex < numCols) {
                    value = 1 - fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
                //rightValue
                columnIndex = scai::common::Math::floor<ValueType>(relativeIndex) + 1;
                if (columnIndex < numCols) {
                    value = fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
            }

            scai::lama::CSRSparseMatrix<ValueType> csrMatrix;
            csrMatrix.allocate(numCols, numColsNew);
            csrMatrix.fillFromAssembly(assembly);

            rMat.swap(csrMatrix);
        }

        /*! \brief Calculates the time step to a corresponding continous time
        \param time continous time in seconds
        \param DT time sampling interval in seconds
        */
        template <typename ValueType>
        scai::IndexType time2index(ValueType time, ValueType DT)
        {
            return (static_cast<scai::IndexType>(time / DT + 0.5));
        }
    }
}
