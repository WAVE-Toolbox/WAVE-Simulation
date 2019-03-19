#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

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

        /*! \brief Calculate the row indices needed for CSR resampling matrix.
        \param csrIA row indices
        \param read_tInt Read Access to query point vector
        \param numCols number of samples before resampling
        \param numColsNew number of samples after resampling
        \param shift number of columns the resampling should be shifted
        */
        template <typename ValueType>
        void calcResampleIA(scai::hmemo::HArray<scai::IndexType> &csrIA, scai::hmemo::ReadAccess<ValueType> &read_tInt, scai::IndexType numCols, scai::IndexType numColsNew, scai::IndexType shift)
        {
            scai::IndexType rawIA[numCols + 1];
            scai::IndexType rawIAProto[numCols + 1]; // number of entries per row

            for (scai::IndexType i = 0; i < numCols + 1; ++i) {
                rawIA[i] = 0;
                rawIAProto[i] = 0;
            }

            scai::IndexType counter;
            for (scai::IndexType i = 0; i < numColsNew; ++i) {
                counter = scai::IndexType(read_tInt[i]) + 1;
                rawIAProto[counter]++;
            }

            rawIA[1 + shift] = rawIAProto[1];

            for (scai::IndexType i = 2 + shift; i < numCols + 1; ++i)
                rawIA[i] = rawIAProto[i - shift] + rawIA[i - 1];

            rawIA[numCols] += numColsNew - rawIA[numCols]; // in some cases add additional ones to last row to prevent size mismatch

            csrIA.setRawData(numCols + 1, rawIA);
        }

        /*! \brief Calculate the column indices needed for CSR resampling matrix.
        \param csrJA column indices
        \param numColsNew number of samples after resampling
        */
        template <typename ValueType>
        void calcResampleJA(scai::hmemo::HArray<scai::IndexType> &csrJA, scai::IndexType numColsNew)
        {
            scai::IndexType rawJA[numColsNew];

            for (scai::IndexType i = 0; i < numColsNew; ++i) {
                rawJA[i] = i;
            }
            csrJA.setRawData(numColsNew, rawJA);
        }

        /*! \brief Calculate a matrix which resamples the columns.
        \param rMat resampling matrix
        \param numCols number of samples in one row
        \param resamplingCoeff resampling coefficient
        \param shift number of columns the resampling should be shifted
        */
        template <typename ValueType>
        void calcResampleMat(scai::lama::CSRSparseMatrix<ValueType> &rMat, scai::IndexType numCols, ValueType resamplingCoeff, scai::IndexType shift)
        {

            SCAI_ASSERT(shift >= 0, "negative shifts are not allowed.")

            scai::IndexType numColsNew = scai::IndexType(scai::common::Math::floor<ValueType>(ValueType(numCols - 1) / ValueType(resamplingCoeff))) + 1; // number of samples after resampling

            scai::lama::DenseVector<ValueType> tInt(scai::lama::linearDenseVector<ValueType>(numColsNew, 0.0, resamplingCoeff));
            scai::hmemo::HArray<ValueType> *tInt_Ptr = &tInt.getLocalValues();
            scai::hmemo::ReadAccess<ValueType> read_tInt(*tInt_Ptr);

            scai::hmemo::HArray<scai::IndexType> csrIA;
            calcResampleIA<ValueType>(csrIA, read_tInt, numCols, numColsNew, shift);
            read_tInt.release();

            scai::hmemo::HArray<scai::IndexType> csrJA;
            calcResampleJA<ValueType>(csrJA, numColsNew);

            scai::hmemo::HArray<ValueType> csrValues(numColsNew, 1.0);

            scai::lama::CSRStorage<ValueType> csrStorage(numCols, numColsNew, csrIA, csrJA, csrValues);
            scai::lama::CSRSparseMatrix<ValueType> csrMatrix(csrStorage);

            rMat.swap(csrMatrix);
        }

        /*! \brief Calculate a vector which contains the query points for interpolation.
        \param resampleVec Result vector
        \param numCols number of samples in one row before interpolation
        \param resamplingCoeff resampling coefficient
        */
        template <typename ValueType>
        void calcResampleVec(scai::lama::DenseVector<ValueType> &resampleVec, scai::IndexType numCols, ValueType resamplingCoeff)
        {
            scai::IndexType numColsNew = scai::IndexType(scai::common::Math::floor<ValueType>(ValueType(numCols - 1) / ValueType(resamplingCoeff))) + 1; // number of samples after resampling
            resampleVec = scai::lama::linearDenseVector<ValueType>(numColsNew, 0.0, resamplingCoeff);
            resampleVec.binaryOpScalar(resampleVec, 1.0, scai::common::BinaryOp::MODULO, false);
            resampleVec.replicate();
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
