#include "Filter.hpp"

template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::init(ValueType dt, scai::IndexType nt) 
{
    zeroPadding = calcZeroPadding(nt);
    scai::IndexType filterLength = nt + zeroPadding;
    transFcn = scai::lama::fill<scai::lama::Vector<ValueType>>(filterLength, 1.0);
    df = 1/(filterLength*dt);
    fNy = 1/(2*dt);
}

template <typename ValueType>
scai::IndexType KITGPI::Filter::Filter<ValueType>::calcZeroPadding(scai::IndexType nt)
{
    ValueType temp = scai::common::Math::log(nt);
    temp /= scai::common::Math::log(2.0);
    temp = scai::common::Math::ceil(temp);
    temp = scai::common::Math::pow(2.0,temp);
    return temp;
}

template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calc(std::string transFcnFmly, std::string filterType, ValueType fc, scai::IndexType order) 
{
    transFcnFmly = std::transform(transFcnFmly.begin(), transFcnFmly.end(), transFcnFmly.begin(), ::tolower);
    if (transFcnFmly == "butterworth") {
        calcButterworthFilt(filterType, fc, order);
    }
    else
    {
        COMMON_THROWEXCEPTION("Invalid transfere function family.");
    }
}

template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthFilt(std::string filterType, ValueType fc, scai::IndexType order)
{
    filterType = std::transform(filterType.begin(), filterType.end(), filterType.begin(), ::tolower);
    
    scai::IndexType nFreq = fNy/df;
    scai::lama::Vector<ValueType> fPos = scai::lama::linearDenseVector<ValueType>(nFreq+1,0.0,df);
    scai::lama::Vector<ValueType> fNeg = scai::lama::linearDenseVector<ValueType>(nFreq-1,-(nFreq-1)*df,-df);
    scai::lama::Vector<ValueType> f = scai::lama::Vector<ValueType>::cat(fPos,fNeg);
    
    if (filterType == "lp") {
        transFcn = f/fc;
    }
    else if (filterType == "hp") {
        transFcn = transFcn.unaryOp(transFcn,scai::common::UnaryOp::RECIPROCAL);
        transFcn *= fc;
    }
    else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
    
    transFcn = scai::lama::pow<ValueType>(transFcn,2.0*order);
    transFcn += 1;
    transFcn.unaryOp(transFcn,scai::common::UnaryOp::RECIPROCAL);
}

template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::Vector<ValueType> &signal)
{
    scai::IndexType len = 2*fNy/df;
    
    scai::lama::Vector<ValueType> padding = scai::lama::fill<scai::lama::Vector<ValueType>>(zeroPadding, 0.0);
    scai::lama::Vector<ValueType> paddedSignal = scai::lama::Vector<ValueType>::cat(signal,padding);
    
    scai::lama::Vector<scai::common::Complex<scai::RealType<ValueType>>> result;
    scai::lama::fft<ValueType>(result,paddedSignal,len);
    
    result *= transFcn;
    
    scai::lama::ifft<ValueType>(result,result,len);
    
    signal = scai::common::Complex<ValueType>::real(result);
    
}


