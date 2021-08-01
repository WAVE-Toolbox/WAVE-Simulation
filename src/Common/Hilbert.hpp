#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
#include <cmath>
#include <scai/lama/fft.hpp>

using namespace scai;
namespace KITGPI
{
    //! \brief Hilbert namespace
    namespace Hilbert
    {
                
        // copy from https://github.com/keziah55/DetectorBank/tree/1772fb22a82f6a3a4f99b992c7a7b16b90cbe01d/src

        /*! Generate an analytic signal where the real part is the input signal and the
        * imaginary part is the Hilbert transform of the input.
        * 
        * Classes which use both FIR and FFT to perform the transform are provided.
        */
        template <typename ValueType>
        class Hilbert 
        {
            public:    
                Hilbert(){};
                ~Hilbert(){};
                
                virtual void hilbert(lama::DenseMatrix<ValueType> &data) = 0;
                virtual void hilbert(lama::DenseVector<ValueType> &data) = 0;
                virtual void calcHilbertCoefficient() = 0;
                void setCoefficientLength(IndexType const setCoefficientLength);
                                
            protected:
                IndexType kernelSize; // Length of coefficient vector
        };
        
        template <typename ValueType>
        class HilbertFFT : public Hilbert<ValueType>
        {
            public:
                HilbertFFT(){};
                ~HilbertFFT(){};
                
                void hilbert(lama::DenseMatrix<ValueType> &data) override;
                void hilbert(lama::DenseVector<ValueType> &data) override;
                void calcHilbertCoefficient() override;
                
            protected:
                typedef common::Complex<RealType<ValueType>> ComplexValueType;
                lama::DenseVector<ComplexValueType> kernel;
                using Hilbert<ValueType>::kernelSize;
        };
              
        template <typename ValueType>
        class HilbertFIR : public Hilbert<ValueType>
        {
            public:
                HilbertFIR(){kernelSize = 10;};
                ~HilbertFIR(){};
                
                void hilbert(lama::DenseMatrix<ValueType> &data) override;
                void hilbert(lama::DenseVector<ValueType> &data) override;
                void calcHilbertCoefficient() override;

            protected:
                
                /*! length of FIR filter (must be odd) */
                IndexType FIRlength;
                /*! kernel size */
                using Hilbert<ValueType>::kernelSize;
                /*! FIR kernel */
                lama::DenseVector<ValueType> kernel;
        };
    }
}
