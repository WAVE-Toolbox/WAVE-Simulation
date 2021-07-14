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
        void searchAndReplace(scai::lama::DenseVector<ValueType> &searchVector, ValueType threshold, ValueType replaceValue, scai::IndexType compareType)
        {
            // needs rework with new LAMA features
            scai::hmemo::HArray<ValueType> *searchVector_Ptr = &searchVector.getLocalValues();
            scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);

            switch (compareType) {
            case 1: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] < threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 2: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] > threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 3: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] <= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 4: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    if (write_searchVector[i] >= threshold) {
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 5: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
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
        void replaceInvalid(scai::lama::DenseVector<ValueType> &searchVector, ValueType replaceValue)
        {
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

        /*! \brief Replaces velocityP (velocityP < VpVsRatio * velocityS) by velocityS
        *
        \param velocityP Input vector
        \param velocityS Input vector
        */
        template <typename ValueType>
        void replaceVpwithVs(scai::lama::DenseVector<ValueType> &velocityP, scai::lama::DenseVector<ValueType> const &velocityS, ValueType VpVsRatio, scai::IndexType compareType)
        {
            scai::hmemo::HArray<ValueType> *searchVector_Ptr = &velocityP.getLocalValues();
            scai::hmemo::WriteAccess<ValueType> write_searchVector(*searchVector_Ptr);
            
            auto read_velocityS = scai::hmemo::hostReadAccess(velocityS.getLocalValues());
            ValueType replaceValue;
            
            switch (compareType) {
            case 1: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
                    replaceValue = VpVsRatio * read_velocityS[i];
                    if (write_searchVector[i] < replaceValue) { 
                        write_searchVector[i] = replaceValue;
                    }
                }
                break;
            }
            case 2: {
                for (scai::IndexType i = 0; i < write_searchVector.size(); ++i) {
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
        scai::IndexType calcNextPowTwo(scai::IndexType nt)
        {
            ValueType temp = scai::common::Math::log(ValueType(nt));
            temp /= scai::common::Math::log(2.0);
            temp = scai::common::Math::ceil(temp);
            temp = scai::common::Math::pow(ValueType(2.0), temp);
            return temp;
        }

        /*! \brief Calculate a matrix which resamples the columns.
        \param rMat resampling matrix
        \param numCols number of samples in one row
        \param resamplingCoeff resampling coefficient
        */
        template <typename ValueType>
        void calcResampleMat(scai::lama::CSRSparseMatrix<ValueType> &rMat, scai::IndexType numCols, ValueType resamplingCoeff)
        {
            scai::lama::MatrixAssembly<ValueType> assembly;

            scai::IndexType numColsNew = scai::IndexType(scai::common::Math::floor<ValueType>(ValueType(numCols - 1) / ValueType(resamplingCoeff))) + 1; // number of samples after resampling

            scai::IndexType columnIndex = 0;
            ValueType value = 0.0;
            ValueType sampleCoeff = 1.0;
            for (scai::IndexType rowIndex = 0; rowIndex < numColsNew; rowIndex++) {

                ValueType relativeIndex = rowIndex * resamplingCoeff;
                //leftValue
                columnIndex = scai::common::Math::floor<ValueType>(relativeIndex);
                if (columnIndex < numCols) {
                    value = 1 - fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
                //rightValue
                columnIndex = scai::common::Math::floor<ValueType>(relativeIndex) + 1;
                if (columnIndex < numCols) {
                    value = fmod(relativeIndex, sampleCoeff);
                    assembly.push(columnIndex, rowIndex, value);
                }
            }

            scai::lama::CSRSparseMatrix<ValueType> csrMatrix;
            csrMatrix.allocate(numCols, numColsNew);
            csrMatrix.fillFromAssembly(assembly);

            rMat.swap(csrMatrix);
        }

        /*! \brief Calculates the time step to a corresponding continous time
        \param time continous time in seconds
        \param DT time sampling interval in seconds
        */
        template <typename ValueType>
        scai::IndexType time2index(ValueType time, ValueType DT)
        {
            return (static_cast<scai::IndexType>(time / DT + 0.5));
        }

        /*! \brief check is seismic or EM
        \param type equationType
        */
        template <typename ValueType>
        bool checkEquationType(std::string type)
        {
            bool isSeismic = false;
            // transform to lower cases
            std::transform(type.begin(), type.end(), type.begin(), ::tolower);
            
            // Assert correctness of input values
            SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");
        
            if (type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0) {
                isSeismic = true;
            } else if (type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0) {
                isSeismic = false;
            }

            return isSeismic;
        }
                
        /*! \brief apply Hilbert transform to data
        *
        \param data Input matrix
        */
        template <typename ValueType>
        void calcHilbert(scai::lama::DenseMatrix<ValueType> &data)
        {
            typedef scai::common::Complex<scai::RealType<ValueType>> ComplexValueType;
            
            scai::IndexType nfLength = Common::calcNextPowTwo<ValueType>(data.getNumColumns());
            auto nsDist = data.getColDistributionPtr();
            auto nfDist = std::make_shared<scai::dmemo::NoDistribution>(nfLength);
            scai::lama::DenseVector<ComplexValueType> h;
            scai::lama::DenseVector<ComplexValueType> fDataTrace;
            scai::lama::DenseMatrix<ComplexValueType> fData;
            
            /* calculation of the vector h */
            h = scai::lama::fill<scai::lama::DenseVector<ComplexValueType>>(nfLength, 0.0);
            if ((2*(nfLength/2))==nfLength) {
                /* nfLength is even */
                h[0]=1.0;
                h[nfLength/2]=1.0;
                for (int i=1;i<(nfLength/2);i++) {
                    h[i] = 2.0;                 
                }
            } else {
                /* nfLength is odd */
                h[0]=1.0;
                for (int i=1;i<((nfLength+1)/2);i++) {
                    h[i]=2.0;                
                }
            }
            
            fData = scai::lama::cast<ComplexValueType>(data);
            fData.resize(data.getRowDistributionPtr(), std::make_shared<scai::dmemo::NoDistribution>(nfLength));

            scai::lama::fft<ComplexValueType>(fData, 1);

            fData.scaleColumns(h);

            fData *= (1.0 / ValueType(nfLength)); // proper fft normalization

            scai::lama::ifft<ComplexValueType>(fData, 1);
            fData.resize(data.getRowDistributionPtr(), data.getColDistributionPtr());
            fData = -fData;
            data = scai::lama::imag(fData);
        }
        
        /*! \brief calculate envelope of the data
        *
        \param data Input matrix
        */
        template <typename ValueType>
        void calcEnvelope(scai::lama::DenseMatrix<ValueType> &data)
        {
            scai::lama::DenseMatrix<ValueType> dataImag = data;
            scai::lama::DenseVector<ValueType> dataTrace;
            calcHilbert(dataImag);
            data.binaryOp(data, scai::common::BinaryOp::MULT, data);
            dataImag.binaryOp(dataImag, scai::common::BinaryOp::MULT, dataImag);
            data.binaryOp(data, scai::common::BinaryOp::ADD, dataImag);
            for (int i=0; i<data.getNumRows(); i++) {
                data.getRow(dataTrace, i);   
                dataTrace = scai::lama::sqrt(dataTrace);
                data.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
            }
        }
        
        /*! \brief calculate instantaneous phase of the data
        *
        \param data Input matrix
        */
        template <typename ValueType>
        void calcInstantaneousPhase(scai::lama::DenseMatrix<ValueType> &data, scai::IndexType const phaseType)
        {
            scai::lama::DenseMatrix<ValueType> dataImag = data;
            scai::lama::DenseVector<ValueType> dataTrace;
            calcHilbert(dataImag);
            dataImag = -dataImag;
            if (phaseType == 1) { // phase wrapped in [-pi/2 pi/2]
                data.binaryOp(dataImag, scai::common::BinaryOp::DIVIDE, data);
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataTrace = scai::lama::atan(dataTrace);
                    data.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
                }
            } else if (phaseType == 2) { // phase wrapped in [-pi pi]
                scai::lama::DenseVector<ValueType> dataImagTrace;
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataImag.getRow(dataImagTrace, i);   
                    for (scai::IndexType tStep = 0; tStep < data.getNumColumns(); tStep++) {
                        ValueType phase = scai::common::Math::atan2(dataImagTrace.getValue(tStep), dataTrace.getValue(tStep)); 
                        dataTrace.setValue(tStep, phase);
                    }
                    data.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
                }
            } else if (phaseType == 3) { // phase unwrapped
                scai::lama::DenseVector<ValueType> dataImagTrace;
                for (int i=0; i<data.getNumRows(); i++) {
                    data.getRow(dataTrace, i);   
                    dataImag.getRow(dataImagTrace, i);   
                    ValueType phase = scai::common::Math::atan2(dataImagTrace.getValue(0), dataTrace.getValue(0)); // phase wrapped in [-pi pi]
                    dataTrace.setValue(0, phase);
                    dataImagTrace.setValue(0, phase);
                    for (scai::IndexType tStep = 1; tStep < data.getNumColumns(); tStep++) {
                        phase = scai::common::Math::atan2(dataImagTrace.getValue(tStep), dataTrace.getValue(tStep)); // phase wrapped in [-pi pi]
                        dataImagTrace.setValue(tStep, phase);
                        ValueType d = phase - dataImagTrace.getValue(tStep-1);
                        d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
                        phase = dataTrace.getValue(tStep-1) + d; // phase unwrapped
                        dataTrace.setValue(tStep, phase);
                    }
                    data.setRow(dataTrace, i, scai::common::BinaryOp::COPY);               
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
        void calcInstantaneousPhaseResidual(scai::lama::DenseMatrix<ValueType> &dataRedsidual, scai::lama::DenseMatrix<ValueType> dataObs, scai::lama::DenseMatrix<ValueType> dataSyn)
        {
            scai::IndexType phaseType = 3;
            Common::calcInstantaneousPhase(dataObs, phaseType);
            Common::calcInstantaneousPhase(dataSyn, phaseType);
            scai::lama::DenseVector<ValueType> dataObsTrace;
            scai::lama::DenseVector<ValueType> dataSynTrace;
            scai::lama::DenseVector<ValueType> dataResTrace1;
            scai::lama::DenseVector<ValueType> dataResTrace2;
            scai::lama::DenseVector<ValueType> dataResTraceMin;
            for (int i=0; i<dataObs.getNumRows(); i++) {
                dataObs.getRow(dataObsTrace, i);   
                dataSyn.getRow(dataSynTrace, i);   
                dataResTrace1 = dataObsTrace - dataSynTrace;
                scai::IndexType phaseShift = 1;
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
                dataRedsidual.setRow(dataResTraceMin, i, scai::common::BinaryOp::COPY); 
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
    }
}
