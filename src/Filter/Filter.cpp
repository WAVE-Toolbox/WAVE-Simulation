#include "Filter.hpp"

/*! \brief Initialize the filter transfere function. In this state nothing is filtered when applying it.
 \param dt Sampling period of the signal the filter should be applied on
 \param nt Length of the signal the filter should be applied on
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::init(ValueType dt, scai::IndexType nt)
{
    zeroPadding = calcZeroPadding(nt);
    scai::IndexType filterLength = zeroPadding+nt;
    transFcn = scai::lama::fill<scai::lama::DenseVector<ValueType>>(filterLength, 1.0);
    df = 1/(filterLength*dt);
    fNyquist = 1/(2*dt);
}

/*! \brief Calculate amount of zeros needed for padding to the next power of two.
 \param nt Length of the signal the filter should be applied on
 */
template <typename ValueType>
scai::IndexType KITGPI::Filter::Filter<ValueType>::calcZeroPadding(scai::IndexType nt)
{
    ValueType temp = scai::common::Math::log(ValueType(nt));
    temp /= scai::common::Math::log(2.0);
    temp = scai::common::Math::ceil(temp);
    temp = scai::common::Math::pow(ValueType(2.0),temp);
    return temp-nt;
}

/*! \brief Calculate the frequency range corresponding to the output of fft.
 \param frequencyVector Vector which should store the frequency range
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcFrequencyVector(scai::lama::DenseVector<ValueType> &frequencyVector)
{
    scai::IndexType nFreq = fNyquist/df;
    scai::lama::DenseVector<ValueType> fPos = scai::lama::linearDenseVector<ValueType>(nFreq+1,0.0,df);
    scai::lama::DenseVector<ValueType> fNeg = scai::lama::linearDenseVector<ValueType>(nFreq-1,-(nFreq-1)*df,-df);
    frequencyVector.cat(fPos,fNeg);
}

/*! \brief Calculate the transfere function of a Butterworth filter.
 \param transFcnFmly Specifies which transfere function type should be used (currently only "butterworth" is possible)
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calc(std::string transFcnFmly, std::string filterType, ValueType fc, scai::IndexType order) 
{
    std::transform(transFcnFmly.begin(), transFcnFmly.end(), transFcnFmly.begin(), ::tolower);
    if (transFcnFmly == "butterworth") {
        calcButterworthFilt(filterType, fc, order);
    }
    else
    {
        COMMON_THROWEXCEPTION("Invalid transfere function family.");
    }
}

/*! \brief Calculate the transfere function of a Butterworth filter.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthFilt(std::string filterType, ValueType fc, scai::IndexType order)
{
    std::transform(filterType.begin(), filterType.end(), filterType.begin(), ::tolower);
    
    scai::lama::DenseVector<ValueType> f;
    calcFrequencyVector(f);
    
    if (filterType == "lp") {
        transFcn = f/fc;
    }
    else if (filterType == "hp") {
        transFcn.unaryOp(transFcn,scai::common::UnaryOp::RECIPROCAL);
        transFcn *= fc;
    }
    else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
    
    transFcn = scai::lama::pow<ValueType>(transFcn,2.0*order);
    transFcn += 1;
    transFcn = scai::lama::sqrt<ValueType>(transFcn);
    transFcn.unaryOp(transFcn,scai::common::UnaryOp::RECIPROCAL);
}

/*! \brief Application of stored filter on the desired signal.
 \param signal Signal that should be filtered
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseVector<ValueType> &signal)
{
    scai::IndexType len = 2*fNyquist/df;
    SCAI_ASSERT_ERROR(signal.size()+zeroPadding == len,"\nFilter is designed for different input length\n\n");
    
    scai::lama::DenseVector<scai::common::Complex<scai::RealType<ValueType>>> signalTemp1;
    scai::lama::DenseVector<scai::common::Complex<scai::RealType<ValueType>>> signalTemp2;
    
    signal.replicate();
    signalTemp1.replicate();
    scai::lama::fft<ValueType>(signalTemp1,signal,len);
    
    scai::lama::DenseVector<scai::common::Complex<scai::RealType<ValueType>>> complexTransfere;
    complexTransfere.buildComplex(transFcn, scai::lama::fill<scai::lama::DenseVector<ValueType>>(transFcn.size(), 0.0));
    signalTemp1 *= complexTransfere;
    signalTemp1 /= (len/2); // proper fft normalization

    signalTemp2.replicate();
    scai::lama::ifft<scai::common::Complex<scai::RealType<ValueType>>>(signalTemp2,signalTemp1,len-zeroPadding);
    
    signal = scai::lama::real(signalTemp2);
}

/*! \brief Application of stored filter on the desired signal.
 \param signal Signal that should be filtered
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseMatrix<ValueType> &signal)
{
        scai::IndexType nTraces = signal.getNumRows();
        scai::lama::DenseVector<ValueType> thisTrace;
        for (scai::IndexType iTrace = 0; iTrace<nTraces; ++iTrace) {
            signal.getRow(thisTrace,iTrace);
            this->apply(thisTrace);
            signal.setRow(thisTrace,iTrace,scai::common::BinaryOp::COPY);
        }
}
// THIS VERSION OF APPLY SHOULD BE USED WHEN LAMA ALLOWS TRUNCATION IN THE FFT FOR MATRICES
// template <typename ValueType>
// void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseMatrix<ValueType> &signal) 
// {
//     scai::IndexType len = 2*fNyquist/df;
//     SCAI_ASSERT_ERROR(signal.getNumColumns()+zeroPadding == len,"\nFilter is designed for different input length\n\n");
//     
//     scai::lama::DenseMatrix<scai::common::Complex<scai::RealType<ValueType>>> signalTemp1;
//     scai::lama::DenseMatrix<scai::common::Complex<scai::RealType<ValueType>>> signalTemp2;
//     
//     scai::dmemo::DistributionPtr no_dist_numTracesGlobal(new scai::dmemo::NoDistribution(signal.getNumColumns()));
//     scai::dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(signal.getNumRows()));
//     signal.redistribute(no_dist_numParameter, no_dist_numTracesGlobal);
// 
//     scai::lama::fft<ValueType>(signalTemp1,signal,0,len);
// 
//     scai::lama::DenseVector<scai::common::Complex<scai::RealType<ValueType>>> complexTransfere;
//     complexTransfere.buildComplex(transFcn, scai::lama::fill<scai::lama::DenseVector<ValueType>>(transFcn.size(), 0.0));
//   
//     signalTemp1.scaleRows(complexTransfere);
//     signalTemp1 *= (2/len); // proper fft normalization
// 
//     scai::lama::ifft<scai::common::Complex<scai::RealType<ValueType>>>(signalTemp2,signalTemp1,0,len-zeroPadding);
//     
//     signal = scai::lama::real(signalTemp2);
// }

template class KITGPI::Filter::Filter<double>;
template class KITGPI::Filter::Filter<float>;
