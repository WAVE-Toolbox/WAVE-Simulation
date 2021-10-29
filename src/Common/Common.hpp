#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
// #include <scai/common/Stencil.hpp>
// #include <scai/lama.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
// #include <scai/lama/matrix/StencilMatrix.hpp>
// #include <scai/lama/matrix/HybridMatrix.hpp>
#include <cmath>
#include <scai/tracing.hpp>
#include <scai/lama/fft.hpp>
#include "../Configuration/Configuration.hpp"
#include "Hilbert.hpp"
#include "medianfilter.hpp"

using namespace scai;
namespace KITGPI
{
    //! \brief Common namespace
    namespace Common
    {
        template <typename T>
        void readColumnFromFile(std::string filename, std::vector<T> &values, unsigned int column)
        {
            std::ifstream is(filename);

            if (is.is_open()) {
                std::string line;
                int lineNumber = 0;
                while (getline(is, line)) {
                    std::vector<int> vec;
                    lineNumber++;
                    std::stringstream strings(line);

                    char firstChar = strings.peek();
                    // check for comment or empty lines or lines with only whitespace or tabs
                    if ((firstChar == '#') || (line.empty() || (std::all_of(line.begin(), line.end(), isspace)))) {
                        continue;
                    } else {
                        T tempStr;
                        while (strings >> tempStr) {
                            vec.push_back(tempStr);
                        }
                        if (vec.size() < column + 1) {
                            COMMON_THROWEXCEPTION("line " << lineNumber << " in file " << filename << " has no column Nr. " << column);
                        }
                        values.push_back(vec[column]);
                    }
                }

            } else {
                COMMON_THROWEXCEPTION("Could not open grid configuration ");
            }
        }

        /*! \brief Searches for all values in searchVector which are related to threshold by relation compareType and replaces them with replaceValue.
        *
        \param searchVector The vector to be searched. 
        \param threshold The threshhold which the values are compared to.
        \param replaceValue The value by which the enties found in searchVector are replaced with.
        \param compareType The relation the values in searchVector and threshold should have. Possible values are 1 := <, 2 := >, 3 := <=, 4 := >=, 5 := ==
        */
        template <typename ValueType>
        void searchAndReplace(lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, IndexType compareType)
        {
            // needs rework with new LAMA features
            hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);

            switch (compareType) {
            case 1: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] < threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 2: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] > threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 3: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] <= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 4: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] >= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 5: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
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
        void replaceInvalid(lama::DenseVector<ValueType> &searchVector, ValueType replaceValue)
        {
            // needs rework with new LAMA features
            hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);

            for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                if (std::isnan(write_searchVector[i]) || write_searchVector[i] == std::numeric_limits<ValueType>::infinity() || -write_searchVector[i] == std::numeric_limits<ValueType>::infinity()) { //std::isinf doesn't work for whatever reason
                    write_searchVector[i] = replaceValue;
                }
            }

            write_searchVector.release();
        }

        /*! \brief Replaces velocityP (velocityP < VpVsRatio * velocityS) by velocityS
        *
        \param velocityP Input vector
        \param velocityS Input vector
        */
        template <typename ValueType>
        void replaceVpwithVs(lama::DenseVector<ValueType> &velocityP, lama::DenseVector<ValueType> const &velocityS, ValueType VpVsRatio, IndexType compareType)
        {
            hmemo::HArray<ValueType> *searchVector_Ptr = &velocityP.getLocalValues();
            hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);
            
            auto read_velocityS = hmemo::hostReadAccess(velocityS.getLocalValues());
            ValueType replaceValue;
            
            switch (compareType) {
            case 1: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    replaceValue = VpVsRatio * read_velocityS[i];
                    if (write_searchVector[i] < replaceValue) { 
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 2: {
                for (IndexType i = 0; i < write_searchVector.size(); ++i) {
                    replaceValue = VpVsRatio * read_velocityS[i];
                    if (write_searchVector[i] > replaceValue) { 
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            default:
                COMMON_THROWEXCEPTION("Invalid compareType. Has to be < 3 but is " << compareType)
            }

            write_searchVector.release();
        }

        /*! \brief Calculate the next power of two.
        \param nt number the next power of two should be calculated for
        */
        template <typename ValueType>
        IndexType calcNextPowTwo(IndexType nt)
        {
            ValueType temp = common::Math::log(ValueType(nt));
            temp /= common::Math::log(2.0);
            temp = common::Math::ceil(temp);
            temp = common::Math::pow(ValueType(2.0), temp);
            return temp;
        }

        /*! \brief Calculate a matrix which resamples the columns.
        \param rMat resampling matrix
        \param numCols number of samples in one row
        \param resamplingCoeff resampling coefficient
        */
        template <typename ValueType>
        void calcResampleMat(lama::CSRSparseMatrix<ValueType> &rMat, IndexType numCols, ValueType resamplingCoeff)
        {
            lama::MatrixAssembly<ValueType> assembly;

            IndexType numColsNew = IndexType(common::Math::floor<ValueType>(ValueType(numCols - 1) / ValueType(resamplingCoeff))) + 1; // number of samples after resampling

            IndexType columnIndex = 0;
            ValueType value = 0.0;
            ValueType sampleCoeff = 1.0;
            for (IndexType rowIndex = 0; rowIndex < numColsNew; rowIndex++) {

                ValueType relativeIndex = rowIndex * resamplingCoeff;
                //leftValue
                columnIndex = common::Math::floor<ValueType>(relativeIndex);
                if (columnIndex < numCols) {
                    value = 1 - fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
                //rightValue
                columnIndex = common::Math::floor<ValueType>(relativeIndex) + 1;
                if (columnIndex < numCols) {
                    value = fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
            }

            lama::CSRSparseMatrix<ValueType> csrMatrix;
            csrMatrix.allocate(numCols, numColsNew);
            csrMatrix.fillFromAssembly(assembly);

            rMat.swap(csrMatrix);
        }

        /*! \brief Calculates the time step to a corresponding continous time
        \param time continous time in seconds
        \param DT time sampling interval in seconds
        */
        template <typename ValueType>
        IndexType time2index(ValueType time, ValueType DT)
        {
            return (static_cast<IndexType>(time / DT + 0.5));
        }

        /*! \brief check is seismic or EM
        \param type equationType
        */
        template <typename ValueType>
        bool checkEquationType(std::string type)
        {
            bool isSeismic = true;
            // transform to lower cases
            std::transform(type.begin(), type.end(), type.begin(), ::tolower);
            
            // Assert correctness of input values
            SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("viscosh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");
        
            if (type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("viscosh") == 0) {
                isSeismic = true;
            } else if (type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0) {
                isSeismic = false;
            }

            return isSeismic;
        }
                        
        /*! \brief calculate envelope of the data
        *
        \param data Input matrix
        */
        template <typename ValueType>
        void calcEnvelope(lama::DenseMatrix<ValueType> &data)
        {
            lama::DenseMatrix<ValueType> dataImag = data;
            lama::DenseVector<ValueType> dataTrace;
            Hilbert::HilbertFFT<ValueType> hilbertHandler;
//             Hilbert::HilbertFIR<ValueType> hilbertHandler;
            IndexType kernelSize = Common::calcNextPowTwo<ValueType>(data.getNumColumns());  
            hilbertHandler.setCoefficientLength(kernelSize);
            hilbertHandler.calcHilbertCoefficient();
            hilbertHandler.hilbert(dataImag);
            data.binaryOp(data, common::BinaryOp::MULT, data);
            dataImag.binaryOp(dataImag, common::BinaryOp::MULT, dataImag);
            data.binaryOp(data, common::BinaryOp::ADD, dataImag);
            for (int i=0; i<data.getNumRows(); i++) {
                data.getRow(dataTrace, i);   
                dataTrace = lama::sqrt(dataTrace);
                data.setRow(dataTrace, i, common::BinaryOp::COPY);               
            }
        }
        
        /*! \brief calculate instantaneous phase of the data
        *
        \param data Input matrix
        */
        template <typename ValueType>
        void calcInstantaneousPhase(lama::DenseMatrix<ValueType> &data, IndexType const phaseType)
        {
            lama::DenseMatrix<ValueType> dataImag = data;
            lama::DenseVector<ValueType> dataTrace;
//             calcHilbertFFT(dataImag);
            dataImag = -dataImag;
            if (phaseType == 1) { // phase wrapped in [-pi/2 pi/2]
                data.binaryOp(dataImag, common::BinaryOp::DIVIDE, data);
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataTrace = lama::atan(dataTrace);
                    data.setRow(dataTrace, i, common::BinaryOp::COPY);               
                }
            } else if (phaseType == 2) { // phase wrapped in [-pi pi]
                lama::DenseVector<ValueType> dataImagTrace;
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataImag.getRow(dataImagTrace, i);   
                    for (IndexType tStep = 0; tStep < data.getNumColumns(); tStep++) {
                        ValueType phase = common::Math::atan2(dataImagTrace.getValue(tStep), dataTrace.getValue(tStep)); 
                        dataTrace.setValue(tStep, phase);
                    }
                    data.setRow(dataTrace, i, common::BinaryOp::COPY);               
                }
            } else if (phaseType == 3) { // phase unwrapped
                lama::DenseVector<ValueType> dataImagTrace;
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataImag.getRow(dataImagTrace, i);   
                    ValueType phase = common::Math::atan2(dataImagTrace.getValue(0), dataTrace.getValue(0)); // phase wrapped in [-pi pi]
                    dataTrace.setValue(0, phase);
                    dataImagTrace.setValue(0, phase);
                    for (IndexType tStep = 1; tStep < data.getNumColumns(); tStep++) {
                        phase = common::Math::atan2(dataImagTrace.getValue(tStep), dataTrace.getValue(tStep)); // phase wrapped in [-pi pi]
                        dataImagTrace.setValue(tStep, phase);
                        ValueType d = phase - dataImagTrace.getValue(tStep-1);
                        d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
                        phase = dataTrace.getValue(tStep-1) + d; // phase unwrapped
                        dataTrace.setValue(tStep, phase);
                    }
                    data.setRow(dataTrace, i, common::BinaryOp::COPY);               
                }
            }
        }
        
        /*! \brief calculate instantaneous phase residual of two data
        *
        \param dataResidual Output data residual
        \param dataObs Input data observed
        \param dataSyn Input data simulated
        */
        template <typename ValueType>
        void calcInstantaneousPhaseResidual(lama::DenseMatrix<ValueType> &dataRedsidual, lama::DenseMatrix<ValueType> dataObs, lama::DenseMatrix<ValueType> dataSyn)
        {
            IndexType phaseType = 3;
            Common::calcInstantaneousPhase(dataObs, phaseType);
            Common::calcInstantaneousPhase(dataSyn, phaseType);
            lama::DenseVector<ValueType> dataObsTrace;
            lama::DenseVector<ValueType> dataSynTrace;
            lama::DenseVector<ValueType> dataResTrace1;
            lama::DenseVector<ValueType> dataResTrace2;
            lama::DenseVector<ValueType> dataResTraceMin;
            for (int i=0; i<dataObs.getNumRows(); i++) {
                dataObs.getRow(dataObsTrace, i);   
                dataSyn.getRow(dataSynTrace, i);   
                dataResTrace1 = dataObsTrace - dataSynTrace;
                IndexType phaseShift = 1;
                dataResTrace2 = dataResTrace1 + phaseShift * 2 * M_PI;
                dataResTraceMin = dataResTrace1;                
                if (dataResTrace2.l2Norm() <= dataResTrace1.l2Norm()) {    
                    while (dataResTrace2.l2Norm() <= dataResTraceMin.l2Norm()) {
                        dataResTraceMin = dataResTrace2;
                        phaseShift += 1;
                        dataResTrace2 = dataResTrace1 + phaseShift * 2 * M_PI;
                    }
                    phaseShift -= 1;
                } else {
                    phaseShift = -1;
                    dataResTrace2 = dataResTrace1 + phaseShift * 2 * M_PI;
                    dataResTraceMin = dataResTrace1;  
                    if (dataResTrace2.l2Norm() <= dataResTrace1.l2Norm()) {               
                        while (dataResTrace2.l2Norm() <= dataResTraceMin.l2Norm()) {
                            dataResTraceMin = dataResTrace2;                
                            phaseShift -= 1;
                            dataResTrace2 = dataResTrace1 + phaseShift * 2 * M_PI;
                        }
                        phaseShift += 1;
                    } else {
                        phaseShift = 0;
                    }
                }
//                 std::cout<< "phaseShift " << i << " = " << phaseShift <<std::endl;
                dataRedsidual.setRow(dataResTraceMin, i, common::BinaryOp::COPY); 
            }
        }
                
        /*! \brief Get the parameter from a stream configuration file
        \param config config handle
        \param parameterName parameter name
        */   
        template <typename ReturnType>
        ReturnType getFromStreamFile(KITGPI::Configuration::Configuration config, std::string const &parameterName)
        {   
            std::string tempName = parameterName;
            ReturnType temp;
            bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
            if (!useStreamConfig) {
                temp = config.get<ReturnType>(tempName);
            } else {
                Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));                
                temp = configBig.get<ReturnType>(tempName);
            }      
            return (temp);
        }        
        
        /*! \brief Apply a 2D median filter to model/gradient parameter to filter out the extreme values
        \param vecter2D model/gradient parameter vector
        \param NX Number of grid points in x-direction
        \param NY Number of grid points in y-direction
        \param spatialFDorder spatial FD order used in Derivatives calculation
        */     
        template<typename ValueType>
        void applyMedianFilterTo2DVector(lama::DenseVector<ValueType> &vecter2D, IndexType NX, IndexType NY, IndexType spatialFDorder)
        {    
            typedef float element;
            element signal_2D[ NY * NX ] = {0};
//             signal_2D = lama::cast<element>(vecter2D);
            for (int i = 0; i < NY * NX; i++) {
                signal_2D[i] = vecter2D[i];
            }
            element *result;
            result = (element *)malloc(NY * NX * sizeof(element));

            medianfilter(signal_2D, result, NX, NY, spatialFDorder*2+1, 0);
            
            for (int i = 0; i < NY * NX; i++) {
                vecter2D[i] = *(result + i);
            }
//             vecter2D = lama::cast<ValueType>(result);
        }
        
        /*! \brief Apply a 1D median filter to model/gradient parameter to filter out the extreme values
        \param vecter1D model/gradient parameter vector
        \param N Number of grid points
        \param spatialFDorder spatial FD order used in Derivatives calculation
        */     
        template<typename ValueType>
        void applyMedianFilterTo1DVector(lama::DenseVector<ValueType> &vecter1D, IndexType spatialFDorder)
        {    
            typedef float element;
            int N = vecter1D.size();
            element signal_1D[N] = {0};
//             signal_1D = lama::cast<element>(vecter1D);
            for (int i = 0; i < N; i++) {
                signal_1D[i] = vecter1D[i];
            }
            element *result;
            result = (element *)malloc(N * sizeof(element));

            medianfilter(signal_1D, result, N, spatialFDorder*2+1);
            
//             vecter1D = lama::cast<ValueType>(result);
            for (int i = 0; i < N; i++) {
                vecter1D[i] = *(result + i);
            }
        }
    }
}
