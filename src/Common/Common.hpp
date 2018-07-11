#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>

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
        template<typename ValueType> 
        void searchAndReplace(scai::lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, scai::IndexType compareType) {
            
            // needs rework with new LAMA features
            scai::hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);
            
            switch(compareType) {
                case 1:
                {
                    for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                        if (write_searchVector[i] < threshold) { 
                            write_searchVector[i] = replaceValue; 
                        }     
                    }
                    break; 
                }
                case 2:
                {
                    for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                        if (write_searchVector[i] > threshold) { 
                            write_searchVector[i] = replaceValue; 
                        }     
                    }
                    break; 
                }
                case 3:
                {
                    for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                        if (write_searchVector[i] <= threshold) { 
                            write_searchVector[i] = replaceValue; 
                        }     
                    }
                    break;
                }
                case 4:
                {
                    for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                        if (write_searchVector[i] >= threshold) { 
                            write_searchVector[i] = replaceValue; 
                        }     
                    }
                    break; 
                }
                case 5:
                {
                    for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                        if (write_searchVector[i] == threshold) { 
                            write_searchVector[i] = replaceValue; 
                        }     
                    }
                    break; 
                }
                default:
                    COMMON_THROWEXCEPTION("Invalid compareType. Has to be < 6 but is "<<compareType)

            }
            
            write_searchVector.release();
        }
        
        /*! \brief Replaces NaN and Inf values in a vector by a given value
        *
        \param searchVector Input vector
        \param replaceValue Value NaN and Inf are set to
        */
        template<typename ValueType> 
        void replaceInvalid(scai::lama::DenseVector<ValueType> &searchVector, ValueType replaceValue) {
            
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

    }
}
