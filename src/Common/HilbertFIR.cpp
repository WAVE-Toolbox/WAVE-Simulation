#include "Hilbert.hpp"

using namespace KITGPI;
using namespace scai;

/*! \brief calculate Hilbert transform coefficient in time domain
*
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFIR<ValueType>::calcHilbertCoefficient()
{
    FIRlength = kernelSize*2-1;
    /* calculation of the vector kernel */
    IndexType m {-kernelSize + 1};
    ValueType k;
    ValueType w;
    ValueType alpha;   
    ValueType a0 = 0.42; // coefficient for blackman function
    ValueType a1 = 0.5;
    ValueType a2 = 0.08;     
    lama::DenseVector<ValueType> temp(kernelSize, 0.0); 
    kernel = temp;
    
    for (IndexType n {0}; n < kernelSize; n++) {
        k = -(cos(m*M_PI)-1) / (m*M_PI);
        alpha = M_PI*n/(kernelSize-1);
        // blackman function
        w = a0 - a1*cos(2.*alpha) + a2 * cos(4.*alpha); 
        
        kernel[n] = k * w;
        kernel[kernelSize-n-1] = -k * w;
        
        m += 2;
    }
}

/*! \brief apply Hilbert transform to wavefield data based on FIR
*
\param data Input wavefield vector
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFIR<ValueType>::hilbert(lama::DenseMatrix<ValueType> &data)
{
    IndexType halfklen {FIRlength / 2};
    IndexType koffset {0};
    IndexType kmin, kmax;
    ValueType h;       
    IndexType signalSize = data.getNumColumns(); 
    lama::DenseVector<ValueType> dataTrace;          
    lama::DenseVector<ValueType> dataTraceHilbert(signalSize, 0);          
    
    if (halfklen == kernelSize) {
        koffset = 1;
    }
    std::cout<< "hilbert matrix kernelSize "<< kernelSize<< std::endl;
    
    for (int itrace=0; itrace<data.getNumRows(); itrace++) {
        data.getRow(dataTrace, itrace);
        for (IndexType n {halfklen}; n < signalSize+halfklen-1; n++) {
            kmin = std::max(static_cast<int>(n-(FIRlength-1)), 0) + koffset;
            kmax = std::min(n, signalSize-1) + koffset;
            
            h = 0.0;
            
            for (IndexType i {kmin}; i < kmax-1; i+=2) {
                h += dataTrace[i] * kernel[(n-i)/2];
            }
            
            dataTraceHilbert[n-halfklen] = h;
        }
        data.setRow(dataTraceHilbert, itrace, scai::common::BinaryOp::COPY);
    }
}

/*! \brief apply Hilbert transform to wavefield data based on FIR
*
\param data Input wavefield vector
*/
template <typename ValueType>
void KITGPI::Hilbert::HilbertFIR<ValueType>::hilbert(lama::DenseVector<ValueType> &data)
{
    IndexType halfklen {FIRlength / 2};
    IndexType koffset {0};
    IndexType kmin, kmax;
    ValueType h;       
    IndexType signalSize = data.size();
    lama::DenseVector<ValueType> dataHilbert = data;          
    
    if (halfklen == kernelSize) {
        koffset = 1;
    }
    std::cout<< "hilbert vector kernelSize "<< kernelSize<< std::endl;
    
    for (IndexType n {halfklen}; n < signalSize+halfklen-1; n++) {
        kmin = std::max(static_cast<int>(n-(FIRlength-1)), 0) + koffset;
        kmax = std::min(n, signalSize-1) + koffset;
        
        h = 0.0;
        
        for (IndexType i {kmin}; i < kmax-1; i+=2) {
            h += data[i] * kernel[(n-i)/2];
        }
        
        dataHilbert[n-halfklen] = h;
    }
    data = dataHilbert;
}

template class KITGPI::Hilbert::HilbertFIR<double>;
template class KITGPI::Hilbert::HilbertFIR<float>;
