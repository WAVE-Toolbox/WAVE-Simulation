#include "Filter.hpp"

/*! \brief Initialize the filter transfere function. In this state nothing is filtered when applying it.
 \param dt Sampling period of the signal the filter should be applied on
 \param nt Length of the signal the filter should be applied on
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::init(ValueType dt, scai::IndexType nt)
{
    SCAI_ASSERT_ERROR(dt != 0.0, "Can't initialize filter with dt = 0.0")
    SCAI_ASSERT_ERROR(nt != 0, "Can't initialize filter with nt = 0")
    zeroPadding = calcZeroPadding(nt);
    scai::IndexType filterLength = zeroPadding + nt;
    transFcn = scai::lama::fill<scai::lama::DenseVector<ValueType>>(filterLength, 1.0);
    df = 1 / (filterLength * dt);
    fNyquist = 1 / (2 * dt);
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
    temp = scai::common::Math::pow(ValueType(2.0), temp);
    return temp - nt;
}

/*! \brief Calculate the frequency range corresponding to the output of fft.
 \param frequencyVector Vector which should store the frequency range
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcFrequencyVector(scai::lama::DenseVector<ValueType> &frequencyVector)
{
    scai::IndexType nFreq = fNyquist / df;
    scai::lama::DenseVector<ValueType> fPos = scai::lama::linearDenseVector<ValueType>(nFreq + 1, 0.0, df);
    scai::lama::DenseVector<ValueType> fNeg = scai::lama::linearDenseVector<ValueType>(nFreq - 1, -(nFreq - 1) * df, df);
    frequencyVector.cat(fPos, fNeg);
}

/*! \brief Calculate the normalized frequency matrix needed for causal Butterworth filters
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc Corner frequency in Hz
 \param frequencyVector Vector which defines the frequency range
 \param frequencyMat Normalized frequency matrix which should be calculated
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcFrequencyMat(std::string filterType, ValueType fc, scai::lama::DenseVector<ValueType> &frequencyVector, scai::lama::DenseMatrix<ComplexValueType> &frequencyMat)
{
    scai::lama::DenseVector<ComplexValueType> fNormCompl;
    scai::lama::DenseVector<ValueType> freqVecNorm;
    
    if (filterType == "lp") {
        freqVecNorm = frequencyVector/fc;
        fNormCompl.buildComplex(scai::lama::fill<scai::lama::DenseVector<ValueType>>(frequencyVector.size(), 0.0), freqVecNorm);
    } else if (filterType == "hp") {
        freqVecNorm = -fc/frequencyVector;
        fNormCompl.buildComplex(scai::lama::fill<scai::lama::DenseVector<ValueType>>(frequencyVector.size(), 0.0), freqVecNorm);
    } else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
    
    for (scai::IndexType iCol = 0; iCol < frequencyMat.getNumColumns(); ++iCol) {
        frequencyMat.setColumn(scai::lama::eval<scai::lama::DenseVector<ComplexValueType>>(scai::lama::pow<ComplexValueType>(fNormCompl,iCol)), iCol, scai::common::BinaryOp::COPY);
    }
    frequencyVector.writeToFile("/home/dkrieger/Masterarbeit/freqVec.mtx");
    frequencyMat.writeToFile("/home/dkrieger/Masterarbeit/freqMat.mtx");
}

/*! \brief Calculate the transfere function of a Butterworth filter.
 \param transFcnFmly Specifies which transfere function type should be used (currently only "butterworth" is possible)
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calc(std::string transFcnFmly, std::string filterType, scai::IndexType order, ValueType fc1, ValueType fc2)
{
    std::transform(transFcnFmly.begin(), transFcnFmly.end(), transFcnFmly.begin(), ::tolower);
    if (transFcnFmly == "butterworth") {
        calcButterworthFilt(filterType, order, fc1, fc2);
    } else {
        COMMON_THROWEXCEPTION("Invalid transfere function family.");
    }
}

/*! \brief Calculate the transfere function of a Butterworth filter.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthFilt(std::string filterType, scai::IndexType order, ValueType fc1, ValueType fc2)
{
    SCAI_ASSERT_ERROR(fc1 > 0, "Lower corner frequency of filter has to be greater than zero.")
    std::transform(filterType.begin(), filterType.end(), filterType.begin(), ::tolower);

    scai::lama::DenseVector<ValueType> freqVec;
    calcFrequencyVector(freqVec);
    
    scai::lama::DenseMatrix<ComplexValueType> freqMat(freqVec.size(), order+1);
    calcFrequencyMat(filterType, fc1, freqVec, freqMat);
    
    if (filterType == "lp") {
        calcButterworthLp(transFcn, freqVec, order, fc1);
    } else if (filterType == "hp") {
        calcButterworthHp(transFcn, freqVec, order, fc1);
    } else if (filterType == "bp") {
        calcButterworthBp(transFcn, freqVec, order, fc1, fc2);
    } else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
}

/*! \brief Calculate the transfere function of a Butterworth low-pass filter.
 \param transFcnTmp Vector where the transfere function gets stored
 \param freqVec Frequency vector
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthLp(scai::lama::DenseVector<ValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc)
{
    transFcnTmp = freqVec / fc;
    transFcnTmp = scai::lama::pow<ValueType>(transFcnTmp, 2.0 * order);
    transFcnTmp += 1;
    transFcnTmp = scai::lama::sqrt<ValueType>(transFcnTmp);
    transFcnTmp.unaryOp(transFcnTmp, scai::common::UnaryOp::RECIPROCAL);
}

/*! \brief Calculate the transfere function of a Butterworth high-pass filter.
 \param transFcnTmp Vector where the transfere function gets stored
 \param freqVec Frequency vector
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthHp(scai::lama::DenseVector<ValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc)
{
    transFcnTmp = fc / freqVec;
    transFcnTmp = scai::lama::pow<ValueType>(transFcnTmp, 2.0 * order);
    transFcnTmp += 1;
    transFcnTmp = scai::lama::sqrt<ValueType>(transFcnTmp);
    transFcnTmp.unaryOp(transFcnTmp, scai::common::UnaryOp::RECIPROCAL);
    transFcnTmp[0] = 0;
}

/*! \brief Calculate the transfere function of a Butterworth band-pass filter.
 \param transFcnTmp Vector where the transfere function gets stored
 \param freqVec Frequency vector
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthBp(scai::lama::DenseVector<ValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc1, ValueType fc2)
{
    SCAI_ASSERT_ERROR(fc2 != 0.0, "Upper corner frequency of band-pass filter can't be zero")
    calcButterworthLp(transFcnTmp, freqVec, order, fc2);
    scai::lama::DenseVector<ValueType> transFcnHp(transFcnTmp);
    calcButterworthHp(transFcnHp, freqVec, order, fc1);
    transFcnTmp *= transFcnHp;
}

/*! \brief Application of stored filter on the desired signal.
 \param signal Signal that should be filtered
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseVector<ValueType> &signal)
{
    scai::IndexType len = 2 * fNyquist / df;
    SCAI_ASSERT_ERROR(signal.size() + zeroPadding == len, "\nFilter is designed for different input length\n\n");

    scai::lama::DenseVector<ComplexValueType> signalTemp1;
    scai::lama::DenseVector<ComplexValueType> signalTemp2;

    signal.replicate();
    signalTemp1.replicate();
    scai::lama::fft<ValueType>(signalTemp1, signal, len);

    scai::lama::DenseVector<ComplexValueType> complexTransfere;
    complexTransfere.buildComplex(transFcn, scai::lama::fill<scai::lama::DenseVector<ValueType>>(transFcn.size(), 0.0));
    signalTemp1 *= complexTransfere;
    signalTemp1 /= (len / 2); // proper fft normalization

    signalTemp2.replicate();
    scai::lama::ifft<ComplexValueType>(signalTemp2, signalTemp1, len - zeroPadding);

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
    for (scai::IndexType iTrace = 0; iTrace < nTraces; ++iTrace) {
        signal.getRow(thisTrace, iTrace);
        this->apply(thisTrace);
        signal.setRow(thisTrace, iTrace, scai::common::BinaryOp::COPY);
    }
}
// THIS VERSION OF APPLY SHOULD BE USED WHEN LAMA ALLOWS TRUNCATION IN THE FFT FOR MATRICES
// template <typename ValueType>
// void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseMatrix<ValueType> &signal)
// {
//     scai::IndexType len = 2*fNyquist/df;
//     SCAI_ASSERT_ERROR(signal.getNumColumns()+zeroPadding == len,"\nFilter is designed for different input length\n\n");
//
//     scai::lama::DenseMatrix<ComplexValueType> signalTemp1;
//     scai::lama::DenseMatrix<ComplexValueType> signalTemp2;
//
//     scai::dmemo::DistributionPtr no_dist_numTracesGlobal(new scai::dmemo::NoDistribution(signal.getNumColumns()));
//     scai::dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(signal.getNumRows()));
//     signal.redistribute(no_dist_numParameter, no_dist_numTracesGlobal);
//
//     scai::lama::fft<ValueType>(signalTemp1,signal,0,len);
//
//     scai::lama::DenseVector<ComplexValueType> complexTransfere;
//     complexTransfere.buildComplex(transFcn, scai::lama::fill<scai::lama::DenseVector<ValueType>>(transFcn.size(), 0.0));
//
//     signalTemp1.scaleRows(complexTransfere);
//     signalTemp1 *= (2/len); // proper fft normalization
//
//     scai::lama::ifft<ComplexValueType>(signalTemp2,signalTemp1,0,len-zeroPadding);
//
//     signal = scai::lama::real(signalTemp2);
// }

template class KITGPI::Filter::Filter<double>;
template class KITGPI::Filter::Filter<float>;
