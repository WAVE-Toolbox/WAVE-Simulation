#include "./Common.hpp"

/*! \brief Searches for all values in searchVector which are related to threshold by relation compareType and replaces them with replaceValue.
 *
 \param searchVector The vector to be searched. 
 \param threshold The threshhold which the values are compared to.
 \param replaceValue The value by which the enties found in searchVector are replaced with.
 \param compareType The relation the values in searchVector and threshold should have. Possible values are 1 := <, 2 := >, 3 := <=, 4 := >=, 5 := ==
 */
template<typename ValueType, typename IndexType> 
void KITGPI::Common::searchAndReplace(scai::lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, IndexType compareType) {
    
    // needs rework with new LAMA features
    scai::hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
    scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);
    
    switch(compareType) {
        case 1:
        {
            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (write_searchVector[i] < threshold) { 
                    write_searchVector[i] = replaceValue; 
                }     
            }
            break; 
        }
        case 2:
        {
            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (write_searchVector[i] > threshold) { 
                    write_searchVector[i] = replaceValue; 
                }     
            }
            break; 
        }
        case 3:
        {
            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (write_searchVector[i] <= threshold) { 
                    write_searchVector[i] = replaceValue; 
                }     
            }
            break;
        }
        case 4:
        {
            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (write_searchVector[i] >= threshold) { 
                    write_searchVector[i] = replaceValue; 
                }     
            }
            break; 
        }
        case 5:
        {
            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
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