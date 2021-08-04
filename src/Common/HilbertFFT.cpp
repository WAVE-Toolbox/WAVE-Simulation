#include "Hilbert.hpp"

using namespace KITGPI;
using namespace scai;

/*! \brief calculate Hilbert transform coefficient in frequency domain
*
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFFT<ValueType>::calcHilbertCoefficient()
{
    /* calculation of the vector kernel */
    kernel = lama::fill<lama::DenseVector<ComplexValueType>>(kernelSize, 0.0);
    if ((2*(kernelSize/2))==kernelSize) {
        /* kernelSize is even */
        kernel[0]=1.0;
        kernel[kernelSize/2]=1.0;
        for (int i=1;i<(kernelSize/2);i++) {
            kernel[i] = 2.0;                 
        }
    } else {
        /* kernelSize is odd */
        kernel[0]=1.0;
        for (int i=1;i<((kernelSize+1)/2);i++) {
            kernel[i]=2.0;                
        }
    }
}

/*! \brief apply Hilbert transform to seismogram data based on FFT
*
\param data Input seismogram matrix
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFFT<ValueType>::hilbert(lama::DenseMatrix<ValueType> &data)
{     
    lama::DenseMatrix<ComplexValueType> fData;
    fData = lama::cast<ComplexValueType>(data);
    fData.resize(data.getRowDistributionPtr(), std::make_shared<dmemo::NoDistribution>(kernelSize));
    lama::fft<ComplexValueType>(fData, 1);
//     std::cout<< "hilbert matrix kernelSize "<< kernelSize<< std::endl;

    fData.scaleColumns(kernel);
    fData *= (1.0 / ValueType(kernelSize)); // proper fft normalization

    lama::ifft<ComplexValueType>(fData, 1);
    fData.resize(data.getRowDistributionPtr(), data.getColDistributionPtr());
    data = lama::imag(fData);
}
        
/*! \brief apply Hilbert transform to wavefield data based on FFT
*
\param data Input wavefield vector
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFFT<ValueType>::hilbert(lama::DenseVector<ValueType> &data)
{
    lama::DenseVector<ComplexValueType> fData;
    fData = lama::cast<ComplexValueType>(data);
    fData.resize(std::make_shared<dmemo::NoDistribution>(kernelSize));
//     std::cout<< "hilbert vector kernelSize "<< kernelSize<< std::endl;
    lama::fft<ComplexValueType>(fData);

    fData *= kernel;
    fData /= kernelSize; // proper fft normalization

    lama::ifft<ComplexValueType>(fData);
    fData.resize(data.getDistributionPtr());
    
    data = lama::imag(fData);
}

template class KITGPI::Hilbert::HilbertFFT<double>;
template class KITGPI::Hilbert::HilbertFFT<float>;
