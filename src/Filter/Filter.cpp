#include "Filter.hpp"

/*! \brief Initialize the filter transfer function. In this state nothing is filtered when applying it.
 \param dt Sampling period of the signal the filter should be applied on
 \param nt Length of the signal the filter should be applied on
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::init(ValueType dt, scai::IndexType nt)
{
    SCAI_ASSERT_ERROR(dt != 0.0, "Can't initialize filter with dt = 0.0")
    SCAI_ASSERT_ERROR(nt != 0, "Can't initialize filter with nt = 0")
    zeroPadding = Common::calcNextPowTwo<ValueType>(nt - 1) - nt;
    scai::IndexType filterLength = zeroPadding + nt;
    filter.buildComplex(scai::lama::fill<scai::lama::DenseVector<ValueType>>(filterLength, 1.0), scai::lama::fill<scai::lama::DenseVector<ValueType>>(filterLength, 0.0));
    df = 1 / (filterLength * dt);
    fNyquist = 1 / (2 * dt);
    NT = nt;
}

/*! \brief Calculate the frequency range corresponding to the output of fft.
 \param frequencyVector Vector which should store the frequency range
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcFrequencyVector(scai::lama::DenseVector<ValueType> &frequencyVector) const
{
    long nFreq = fNyquist / df;
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
void KITGPI::Filter::Filter<ValueType>::calcFrequencyMat(std::string filterType, scai::IndexType order, ValueType fc, scai::lama::DenseVector<ValueType> &frequencyVector, scai::lama::DenseMatrix<ComplexValueType> &frequencyMat)
{
    frequencyMat.allocate(frequencyVector.size(), order + 1);

    scai::lama::DenseVector<ComplexValueType> fNormCompl;
    scai::lama::DenseVector<ValueType> freqVecNorm;

    if (filterType == "lp") {
        freqVecNorm = frequencyVector / fc;
        fNormCompl.buildComplex(scai::lama::fill<scai::lama::DenseVector<ValueType>>(frequencyVector.size(), 0.0), freqVecNorm);
    } else if (filterType == "hp") {
        freqVecNorm = -fc / frequencyVector;
        fNormCompl.buildComplex(scai::lama::fill<scai::lama::DenseVector<ValueType>>(frequencyVector.size(), 0.0), freqVecNorm);
    } else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }

    for (scai::IndexType iCol = 0; iCol <= order; ++iCol) {
        frequencyMat.setColumn(scai::lama::eval<scai::lama::DenseVector<ComplexValueType>>(scai::lama::pow<ComplexValueType>(fNormCompl, iCol)), iCol, scai::common::BinaryOp::COPY);
    }
}

/*! \brief Calculate the coefficients of the Butterworth polynomial in case the order is even
 \param order Filter order
 \param poly Result where the coefficients are stored
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterPoly(scai::IndexType order, scai::lama::DenseVector<ValueType> &poly)
{
    // see https://www.cnblogs.com/xpvincent/p/5557659.html for details (in Chinese)
    poly = scai::lama::fill<scai::lama::DenseVector<ValueType>>(order + 1, 0.0);
    poly[0] = 1.0; //using this as initial factor for the convolution gives the right order of coefficients

    scai::IndexType fPolyLen = Common::calcNextPowTwo<ValueType>(order + 1);

    auto fPolyDist = std::make_shared<scai::dmemo::NoDistribution>(fPolyLen);
    auto polyDist = poly.getDistributionPtr();

    scai::lama::DenseVector<ComplexValueType> fPoly;
    fPoly = scai::lama::cast<ComplexValueType>(poly);
    fPoly.resize(fPolyDist);
    scai::lama::fft<ComplexValueType>(fPoly);

    scai::lama::DenseVector<ValueType> polyTemp = scai::lama::fill<scai::lama::DenseVector<ValueType>>(fPolyLen, 0.0); // Single factor of Butterworth polynomial
    polyTemp[0] = 1.0;
    polyTemp[2] = 1.0;

    scai::lama::DenseVector<ComplexValueType> fPolyTemp;

    for (scai::IndexType iOrder = 1; iOrder <= order / 2; ++iOrder) {
        polyTemp[1] = -2.0 * scai::common::Math::cos((2.0 * iOrder + order - 1.0) / (2.0 * order) * M_PI);
        fPolyTemp = scai::lama::cast<ComplexValueType>(polyTemp);
        fPolyTemp.resize(fPolyDist);
        scai::lama::fft<ComplexValueType>(fPolyTemp);
        fPoly *= fPolyTemp; // convolution of polynomial factors is multiplication in frequency domain
    }

    if (order % 2 != 0) {
        polyTemp[1] = 1.0;
        polyTemp[2] = 0.0;
        fPolyTemp = scai::lama::cast<ComplexValueType>(polyTemp);
        fPolyTemp.resize(fPolyDist);
        scai::lama::fft<ComplexValueType>(fPolyTemp);
        fPoly *= fPolyTemp;
    }

    fPoly /= fPolyLen;
    scai::lama::ifft<ComplexValueType>(fPoly);
    fPoly.resize(polyDist);
    poly = scai::lama::real(fPoly);
}

/*! \brief Calculate the transfer function of a filter.
 \param transFcnFmly Specifies which transfer function type should be used.
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
    } else if (transFcnFmly == "ideal") {
        calcIdealFilt(filterType, order, fc1);
    } else {
        COMMON_THROWEXCEPTION("Invalid transfer function family.");
    }
}

/*! \brief Calculate the transfer function of a Butterworth filter.
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

    if (filterType == "lp") {
        calcButterworthLp(filter, freqVec, order, fc1);
    } else if (filterType == "hp") {
        calcButterworthHp(filter, freqVec, order, fc1);
    } else if (filterType == "bp") {
        calcButterworthBp(filter, freqVec, order, fc1, fc2);
    } else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
}

/*! \brief Calculate the transfer function of a Butterworth low-pass filter.
 \param transFcnTmp Vector where the transfer function gets stored
 \param freqVec Frequency vector
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthLp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc)
{
    scai::lama::DenseMatrix<ComplexValueType> freqMat;
    calcFrequencyMat("lp", order, fc, freqVec, freqMat);

    scai::lama::DenseVector<ValueType> poly;
    calcButterPoly(order, poly);

    scai::lama::DenseVector<ComplexValueType> cPoly;
    cPoly = scai::lama::cast<ComplexValueType>(poly);

    transFcnTmp = freqMat * cPoly;
    transFcnTmp.unaryOp(transFcnTmp, scai::common::UnaryOp::RECIPROCAL);
    transFcnTmp.setValue(0, ComplexValueType(1.0, 0.0));
}

/*! \brief Calculate the transfer function of a Butterworth high-pass filter.
 \param transFcnTmp Vector where the transfer function gets stored
 \param freqVec Frequency vector
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthHp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc)
{
    scai::lama::DenseMatrix<ComplexValueType> freqMat;
    calcFrequencyMat("hp", order, fc, freqVec, freqMat);

    scai::lama::DenseVector<ValueType> poly;
    calcButterPoly(order, poly);

    scai::lama::DenseVector<ComplexValueType> cPoly;
    cPoly = scai::lama::cast<ComplexValueType>(poly);

    transFcnTmp = freqMat * cPoly;
    transFcnTmp.unaryOp(transFcnTmp, scai::common::UnaryOp::RECIPROCAL);
    transFcnTmp.setValue(0, ComplexValueType(0.0, 0.0));
}

/*! \brief Calculate the transfer function of a Butterworth band-pass filter.
 \param transFcnTmp Vector where the transfer function gets stored
 \param freqVec Frequency vector
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcButterworthBp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::lama::DenseVector<ValueType> &freqVec, scai::IndexType order, ValueType fc1, ValueType fc2)
{
    SCAI_ASSERT_ERROR(fc2 != 0.0, "Upper corner frequency of band-pass filter can't be zero")
    calcButterworthLp(transFcnTmp, freqVec, order, fc2);
    scai::lama::DenseVector<ComplexValueType> transFcnHp(transFcnTmp);
    calcButterworthHp(transFcnHp, freqVec, order, fc1);
    transFcnTmp *= transFcnHp;
}

/*! \brief Calculate the transfer function of an ideal filter.
 \param filterType Type of filter: "lp" = low pass, "hp" = high pass
 \param fc1 Lower corner frequency in Hz
 \param fc2 Upper corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcIdealFilt(std::string filterType, scai::IndexType order, ValueType fc)
{
    SCAI_ASSERT_ERROR(fc > 0, "Lower corner frequency of filter has to be greater than zero.")
    std::transform(filterType.begin(), filterType.end(), filterType.begin(), ::tolower);

    if (filterType == "bp") {
        calcIdealBp(filter, order, fc);
    } else {
        COMMON_THROWEXCEPTION("Invalid filter type.");
    }
}

/*! \brief Calculate the transfer function of an ideal band-pass filter (one frequency).
 \param transFcnTmp Vector where the transfer function gets stored
 \param freqVec Frequency vector
 \param fc Corner frequency in Hz
 \param order Filter order
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::calcIdealBp(scai::lama::DenseVector<ComplexValueType> &transFcnTmp, scai::IndexType order, ValueType fc)
{
    IndexType fc1Ind = ceil(fc / df);
    if (order == 0) { // FFT
        scai::IndexType len = 2 * fNyquist / df;
        IndexType fc2Ind = len - fc1Ind;
        // transFcnTmp (filter) has been initialized to 1+0i
        transFcnTmp.setValue(fc1Ind, ComplexValueType(0.0, 0.0));
        transFcnTmp.setValue(fc2Ind, ComplexValueType(0.0, 0.0));
        transFcnTmp = ComplexValueType(1.0, 0.0) - transFcnTmp;
    } else { // DFT
        scai::dmemo::DistributionPtr distNT = std::make_shared<scai::dmemo::NoDistribution>(NT);
        scai::dmemo::DistributionPtr dist1 = std::make_shared<scai::dmemo::NoDistribution>(1);
        L.allocate(distNT, dist1);
        Linv.allocate(dist1, distNT);
        ValueType dt = 1 / (2 * fNyquist);
        ComplexValueType j(0.0, 1.0); 
        scai::lama::DenseVector<ValueType> t = scai::lama::linearDenseVector<ValueType>(NT, 0.0, dt);
        scai::lama::DenseVector<ComplexValueType> temp;
        temp = scai::lama::cast<ComplexValueType>(t);
        temp *= -j * 2.0 * M_PI * fc1Ind * df;
        temp.unaryOp(temp, common::UnaryOp::EXP);
        L.setColumn(temp, 0, common::BinaryOp::COPY);
        temp = scai::lama::cast<ComplexValueType>(t);
        temp *= j * 2.0 * M_PI * fc1Ind * df;
        temp.unaryOp(temp, common::UnaryOp::EXP);
        Linv.setRow(temp, 0, common::BinaryOp::COPY);
    }
}

/*! \brief Application of stored filter on the desired signal.
 \param signal Signal that should be filtered
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseVector<ValueType> &signal) const
{
    scai::IndexType len = 2 * fNyquist / df;
    SCAI_ASSERT_ERROR(signal.size() + zeroPadding == len, "\nFilter is designed for different input length\n\n");

    scai::lama::DenseVector<ComplexValueType> fSignal;
    fSignal = scai::lama::cast<ComplexValueType>(signal);

    if (L.getNumColumns() == 0) {
        fSignal.resize(std::make_shared<scai::dmemo::NoDistribution>(len));
        scai::lama::fft<ComplexValueType>(fSignal);

        fSignal *= filter;
        fSignal /= len; // proper fft normalization

        scai::lama::ifft<ComplexValueType>(fSignal);
        fSignal.resize(signal.getDistributionPtr());
    } else {
        scai::lama::DenseVector<ComplexValueType> temp;
        scai::lama::DenseMatrix<ComplexValueType> Ltemp;
        Ltemp = transpose(L);
        temp = Ltemp * fSignal;
        Ltemp = transpose(Linv);
        fSignal = Ltemp * temp;
        fSignal *= (1.0 / ValueType(NT)); // proper fft normalization
    }
    signal = scai::lama::real(fSignal);
}

/*! \brief Application of stored filter on the desired signal.
 \param signal Signal that should be filtered
 */
template <typename ValueType>
void KITGPI::Filter::Filter<ValueType>::apply(scai::lama::DenseMatrix<ValueType> &signal) const
{
    scai::IndexType len = 2 * fNyquist / df;
    SCAI_ASSERT_ERROR(signal.getNumColumns() + zeroPadding == len, "\nFilter is designed for different input length\n\n");

    scai::lama::DenseMatrix<ComplexValueType> fSignal;
    fSignal = scai::lama::cast<ComplexValueType>(signal);
    
    if (L.getNumColumns() == 0) {
        fSignal.resize(signal.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(len));

        scai::lama::fft<ComplexValueType>(fSignal, 1);

        fSignal.scaleColumns(filter);

        fSignal *= (1.0 / ValueType(len)); // proper fft normalization

        scai::lama::ifft<ComplexValueType>(fSignal, 1);
        fSignal.resize(signal.getRowDistributionPtr(), signal.getColDistributionPtr());
    } else {
        scai::lama::DenseMatrix<ComplexValueType> temp;
        temp = fSignal * L;
        temp *= (2.0 / ValueType(len)); // proper dft normalization for single frequency
        fSignal = temp * Linv;
    }
    signal = scai::lama::real(fSignal);
    if (L.getNumColumns() != 0) {
        signal *= 1.0 / signal.maxNorm(); // scale to maximum amplitude
    }
}

template class KITGPI::Filter::Filter<double>;
template class KITGPI::Filter::Filter<float>;
