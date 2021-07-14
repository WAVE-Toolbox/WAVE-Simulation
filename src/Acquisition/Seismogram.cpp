#include "Seismogram.hpp"
#include "../IO/IO.hpp"
#include "../IO/SUIO.hpp"

#include <scai/utilskernel/HArrayUtils.hpp>
using namespace scai;

//! \brief copy constructor
/*!
 *
 \param rhs seismogram to copy
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>::Seismogram(const Seismogram &rhs)
{
    outputDT = rhs.outputDT;
    DT = rhs.DT;
    seismoType = rhs.seismoType;
    coordinates1D = rhs.coordinates1D;
    sourceIndex = rhs.sourceIndex;
    data = rhs.data;
    resampleMat = rhs.resampleMat;
    outputInstantaneous = rhs.outputInstantaneous;
    frequencyAGC = rhs.frequencyAGC;
    inverseAGC = rhs.inverseAGC;
}
//! \brief swap function
/*!
 *
 * swaps all members from rhs and lhs
 \param rhs seismogram to swap with
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::swap(KITGPI::Acquisition::Seismogram<ValueType> &rhs)
{
    std::swap(DT, rhs.DT);
    std::swap(outputDT, rhs.outputDT);
    std::swap(seismoType, rhs.seismoType);
    std::swap(coordinates1D, rhs.coordinates1D);
    std::swap(sourceIndex, rhs.sourceIndex);
    std::swap(outputInstantaneous, rhs.outputInstantaneous);
    std::swap(frequencyAGC, rhs.frequencyAGC);
    data.swap(rhs.data);
    inverseAGC.swap(rhs.inverseAGC);
}

//! \brief Setter method for the context ptr
/*!
 *
 * This method sets the Context to the coordinates1D and to the seismogram data.
 \param ctx Set ContextPtr
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setContextPtr(scai::hmemo::ContextPtr ctx)
{
    data.setContextPtr(ctx);
    coordinates1D.setContextPtr(ctx);
    inverseAGC.setContextPtr(ctx);
}

//! \brief  writes the seismogram to disk
/*!
 *
 * This method writes the seismogram data to disk. \n
 * The output format is determined by the input parameter `SeismogramFormat` \n
 \param filename base filename of the seismogram
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param seismogramFormat =1 MTX: MatrixMaker format, =4 SU: SeismicUnix format
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::write(scai::IndexType const seismogramFormat, std::string const &filename, Coordinates<ValueType> const &modelCoordinates) const
{
    if (data.getNumValues() > 0) {
        scai::lama::DenseMatrix<ValueType> dataResample;
        dataResample = data * resampleMat;
        
        if (outputInstantaneous == 1) {
            Common::calcEnvelope(dataResample);
        } else if (outputInstantaneous == 2) {
            scai::IndexType phaseType = 2;
            Common::calcInstantaneousPhase(dataResample, phaseType);
        } 

        scai::IndexType seismoFormat = seismogramFormat;
        std::string filenameTmp = filename + "." + SeismogramTypeString[getTraceType()];
        if (seismogramFormat == 5) {
            seismoFormat = 1;
            dataResample = inverseAGC;
            filenameTmp += ".inverseAGC";
        }
        
        switch (seismoFormat) {
        case 4:
            SUIO::writeSU(filenameTmp, dataResample, coordinates1D, outputDT, sourceIndex, modelCoordinates);
            break;
        default:
            IO::writeMatrix(dataResample, filenameTmp, seismoFormat);
            break;
        }
    }
}

//! \brief Read the seismogram from disk
/*!
 *
 * This method reads the seismogram data from disk. \n
\param seismogramFormat =1 MTX: MatrixMaker format, =4 SU: SeismicUnix format
 \param filename Filename to read seismogram
 \param copyDist Boolean: 0 = read data undistributed (default), data is replicated on each process // 1 = read data with existing distribution of data
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::read(scai::IndexType const seismogramFormat, std::string const &filename, bool copyDist)
{
    std::string filenameTmp = filename + "." + SeismogramTypeString[getTraceType()];
    scai::IndexType seismoFormat = seismogramFormat;
    scai::lama::DenseMatrix<ValueType> dataTemp(data);
    if (seismogramFormat == 5) {
        seismoFormat = 1;
        filenameTmp += ".inverseAGC";
        useAGC = true;
    }
    
    switch (seismoFormat) {
    case 4:
        SUIO::readDataSU(filenameTmp, dataTemp, this->getNumSamples(), this->getNumTracesGlobal());
        break;
    default:
        IO::readMatrix(dataTemp, filenameTmp, seismoFormat);
        break;
    }
    if (seismogramFormat == 5) {
        inverseAGC = dataTemp;
    } else {
        data = dataTemp;
    }
}

//! \brief Normalize the seismogram-traces
/*!
 *
 * This method normalized the traces of the seismogram after the time stepping.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::normalizeTrace(scai::IndexType normalizeTraces)
{    
    if (data.getNumValues() > 0 && normalizeTraces > 0) {
        scai::hmemo::HArray<ValueType> tempRow;
        if (normalizeTraces == 1) {
            for (IndexType i = 0; i < getNumTracesLocal(); i++) {
                data.getLocalStorage().getRow(tempRow, i);
                ValueType tempMax = scai::utilskernel::HArrayUtils::maxNorm(tempRow);
                if (tempMax == 0) tempMax = 1;
                scai::utilskernel::HArrayUtils::setScalar(tempRow, tempMax, scai::common::BinaryOp::DIVIDE);
                data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
            }
        } else if (normalizeTraces == 2) {
            for (IndexType i = 0; i < getNumTracesLocal(); i++) {
                data.getLocalStorage().getRow(tempRow, i);
                ValueType tempMax = scai::utilskernel::HArrayUtils::l2Norm(tempRow);
                if (tempMax == 0) tempMax = 1;
                scai::utilskernel::HArrayUtils::setScalar(tempRow, tempMax, scai::common::BinaryOp::DIVIDE);
                data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
            }
        } else if (normalizeTraces == 3 && useAGC) { // normalized by AGC
            data.binaryOp(data, common::BinaryOp::MULT, inverseAGC);
            // normalized to the l2 norm
            for (IndexType i = 0; i < getNumTracesLocal(); i++) {
                data.getLocalStorage().getRow(tempRow, i);
                ValueType tempMax = scai::utilskernel::HArrayUtils::l2Norm(tempRow);
                if (tempMax == 0) tempMax = 1;
                scai::utilskernel::HArrayUtils::setScalar(tempRow, tempMax, scai::common::BinaryOp::DIVIDE);
                data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
            }
            useAGC = false;
        }
    }
}

//! \brief Calculate and get the AGC sum function
/*!
 *
 * This method calculate and get the AGC sum function.
 */
template <typename ValueType>
scai::lama::DenseMatrix<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::getAGCSum()
{    
    scai::lama::DenseMatrix<ValueType> AGCSum;
    if (data.getNumValues() > 0) {
        AGCSum = data;
        AGCSum.scale(0);
        scai::IndexType NAGC = round(2.0 / (frequencyAGC * DT));
        scai::hmemo::HArray<ValueType> tempRow;
        for (IndexType i = 0; i < getNumTracesLocal(); i++) {
            data.getLocalStorage().getRow(tempRow, i);
            scai::hmemo::HArray<ValueType> AGCSumRow(tempRow);
            ValueType var = 0;
            ValueType sumTemp = 0;
            ValueType NWIN = NAGC;
            for (IndexType tStep = getNumSamples()-NAGC; tStep < getNumSamples(); tStep++) {
                var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep);
                sumTemp += var;
            }
            for (IndexType tStep = getNumSamples()-1; tStep >= 0; tStep--) {
                // we have to use decreasing time step because increasing time step generates many negative values in sumTemp which may be caused by zero parts in simulated data.
                if (tStep >= getNumSamples()-NAGC) { // ramping on
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep-NAGC);
                    sumTemp += var;
                    NWIN += 1;
                } else if (tStep >= NAGC && tStep < getNumSamples()-NAGC) { // middle range -- full rms window 
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep-NAGC);
                    sumTemp += var;
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep+NAGC);
                    sumTemp -= var;
                } else if  (tStep < NAGC) { // ramping off 
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep+NAGC);
                    sumTemp -= var;
                    NWIN -= 1;
                }
                ValueType rmsTemp = sumTemp / NWIN;
                scai::utilskernel::HArrayUtils::setVal(AGCSumRow, tStep, rmsTemp, scai::common::BinaryOp::COPY);
            }
            AGCSum.getLocalStorage().setRow(AGCSumRow, i, scai::common::BinaryOp::COPY);
        }        
    }
    return AGCSum;
}
//! \brief Calculate the inverse of AGC function
/*!
 *
 * This method calculate the inverse of AGC function.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::calcInverseAGC()
{    
    if (data.getNumValues() > 0) {
        useAGC = true;
        inverseAGC = data;
        inverseAGC.scale(0);
        scai::IndexType NAGC = round(2.0 / (frequencyAGC * DT));
        scai::hmemo::HArray<ValueType> tempRow;
        for (IndexType i = 0; i < getNumTracesLocal(); i++) {
            data.getLocalStorage().getRow(tempRow, i);
            scai::hmemo::HArray<ValueType> inverseAGCRow(tempRow);
            ValueType var = 0;
            ValueType sumTemp = 0;
            ValueType NWIN = NAGC;
            ValueType waterLevel = scai::utilskernel::HArrayUtils::l2Norm(tempRow);
            waterLevel *= waterLevel;
            waterLevel /= getNumSamples();
            if (waterLevel != 0) {
                waterLevel *= 1e-3;
            } else {
                waterLevel = 1;
            }
            for (IndexType tStep = getNumSamples()-NAGC; tStep < getNumSamples(); tStep++) {
                var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep);
                sumTemp += var * var;
                sumTemp += waterLevel;
            }
            for (IndexType tStep = getNumSamples()-1; tStep >= 0; tStep--) {
                // we have to use decreasing time step because increasing time step generates many negative values in sumTemp which may be caused by zero parts in simulated data.
                if (tStep >= getNumSamples()-NAGC) { // ramping on
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep-NAGC);
                    sumTemp += var * var;
                    sumTemp += waterLevel;
                    NWIN += 1;
                } else if (tStep >= NAGC && tStep < getNumSamples()-NAGC) { // middle range -- full rms window 
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep-NAGC);
                    sumTemp += var * var;
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep+NAGC);
                    sumTemp -= var * var;
                } else if  (tStep < NAGC) { // ramping off 
                    var = scai::utilskernel::HArrayUtils::getVal(tempRow, tStep+NAGC);
                    sumTemp -= var * var;
                    sumTemp -= waterLevel;
                    NWIN -= 1;
                }
                ValueType rmsTemp = sumTemp / NWIN;
                rmsTemp = sqrt(rmsTemp);
                rmsTemp = 1 / rmsTemp;
                scai::utilskernel::HArrayUtils::setVal(inverseAGCRow, tStep, rmsTemp, scai::common::BinaryOp::COPY);
            }
            inverseAGC.getLocalStorage().setRow(inverseAGCRow, i, scai::common::BinaryOp::COPY);
        }        
    }
}

//! \brief Get the inverse of AGC function
/*!
 *
 * This method get the inverse of AGC function.
 */
template <typename ValueType>
scai::lama::DenseMatrix<ValueType> const &KITGPI::Acquisition::Seismogram<ValueType>::getInverseAGC() const
{    
    return inverseAGC;
}

//! \brief Set the inverse of AGC function
/*!
 *
 * This method set the inverse of AGC function.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setInverseAGC(scai::lama::DenseMatrix<ValueType> setInverseAGC)
{    
    inverseAGC = setInverseAGC;
    useAGC = true;
}

//! \brief Normalize the seismogram-traces
/*!
 *
 * This method normalized the traces of the seismogram after the time stepping.
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::getTraceL2norm()
{    
    scai::lama::DenseVector<ValueType> tempRow;
    scai::lama::DenseVector<ValueType> traceL2Norm;
    traceL2Norm.allocate(data.getRowDistributionPtr());
    traceL2Norm = 1.0; // in this state the taper does nothing when applied
    traceL2Norm.setContextPtr(data.getContextPtr());
    
    if (data.getNumValues() > 0) {
        for (IndexType i = 0; i < getNumTracesGlobal(); i++) {
            data.getRow(tempRow, i);
            traceL2Norm.setValue(i, tempRow.l2Norm());
            if (traceL2Norm.getValue(i)==0) traceL2Norm.setValue(i,1);
        }
    }
    return traceL2Norm;
}

//! \brief Normalize the seismogram-traces
/*!
 *
 * This method normalized the traces of the seismogram after the time stepping.
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::getTraceSum()
{    
    scai::lama::DenseVector<ValueType> tempRow;
    scai::lama::DenseVector<ValueType> traceSum;
    traceSum.allocate(data.getRowDistributionPtr());
    traceSum = 1.0; // in this state the taper does nothing when applied
    traceSum.setContextPtr(data.getContextPtr());
    
    if (data.getNumValues() > 0) {
        for (IndexType i = 0; i < getNumTracesGlobal(); i++) {
            data.getRow(tempRow, i);
            traceSum.setValue(i, tempRow.sum());
            if (traceSum.getValue(i)==0) traceSum.setValue(i,1);
        }
    }
    return traceSum;
}

//! \brief Integrate the seismogram-traces
/*!
 *
 * This method integrate the traces of the seismogram.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::integrateTraces()
{
    if (data.getNumValues() > 0) {
        scai::hmemo::HArray<ValueType> tempRow;
        for (IndexType i = 0; i < getNumTracesLocal(); i++) {
            data.getLocalStorage().getRow(tempRow, i);
            for (IndexType j = 0; j < tempRow.size() - 1; j++) {
                tempRow[j + 1] = tempRow[j + 1] * DT + tempRow[j];
            }
            data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
        }
    }
}

//! \brief Filter the seismogram-traces
/*!
 *
 * This method filters the traces of the seismogram.
 \param freqFilter filter object
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::filterTraces(Filter::Filter<ValueType> const &freqFilter)
{
    if (this->getNumSamples() != 0) {
        freqFilter.apply(data);
    }
}

//! \brief Check seismogram for infinite/NaN value
/*!
 *
 * Check last sample of each trace of seismogram for infinite/NaN value
 */
template <typename ValueType>
bool KITGPI::Acquisition::Seismogram<ValueType>::isFinite()
{
    bool result_isfinite=true;
    for (IndexType loc_vals=this->getNumSamples()-1;loc_vals<data.getLocalStorage().getData().size()-1;loc_vals=loc_vals+this->getNumSamples()) {
        if (isfinite(data.getLocalStorage().getData()[loc_vals])==false){
            result_isfinite=false;
            break;
        }
    }
    return(result_isfinite);
}

//! \brief Setter method for the temporal sampling DT
/*!
 *
 * This method will set the temporal sampling DT to this class.
 \param newDT Temporal sampling which will be set in seconds
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setDT(ValueType newDT)
{
    SCAI_ASSERT(newDT >= 0, " DT is smaller zero. ");
    DT = newDT;
}

//! \brief Setter method for the source Index
/*!
 *
 * This method sets the source Index to this class. The source Index will be used for calculation of the header information (e.g. Offset), mainly during seismogram output to disk.
 \param sourceCoord Source Index
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setSourceCoordinate(IndexType sourceIdx)
{
    SCAI_ASSERT_DEBUG(sourceIdx >= 0, "sourceCoord is not valid");
    sourceIndex = sourceIdx;
}

//! \brief Setter method for the #SeismogramType
/*!
 *
 * This method will set the #SeismogramType to this class.
 \param trace Trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setTraceType(SeismogramType trace)
{
    seismoType = trace;
};

//! \brief Setter function for the 1D coordinates of the traces
/*!
 *
 * This method will set the coordinates1D of the traces to this class. 
 * The size of the Index vector has to be equal to the number of global traces.
 \param indeces DenseVector with coordinates1D
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setCoordinates(scai::lama::DenseVector<IndexType> const &indeces)
{
    SCAI_ASSERT_ERROR(indeces.size() == getNumTracesGlobal(), "Given traceType vector has wrong format");
    coordinates1D = indeces;
};

//! \brief Setter method to set matrix for resampling this seismogram.
/*!
 \param rMat Resampling matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setSeismoDT(ValueType seismoDT)
{
    outputDT = seismoDT;
    ValueType resampleCoeff = seismoDT / DT;

    if (this->getNumSamples() != 0) {
        Common::calcResampleMat<ValueType>(resampleMat, getNumSamples(), resampleCoeff);
    }
}

//! \brief Setter method to set outputInstantaneous.
/*!
 \param outputInstantaneous outputInstantaneous
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setInstantaneousTrace(scai::IndexType instantaneousTraces)
{
    outputInstantaneous = instantaneousTraces;
}

//! \brief Setter method to set frequencyAGC.
/*!
 \param setFrequencyAGC setFrequencyAGC
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setFrequencyAGC(ValueType setFrequencyAGC)
{
    frequencyAGC = setFrequencyAGC;
}

//! \brief Getter method for #SeismogramType
/*!
 *
 * This method returns the #SeismogramType of this seismogram.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramType KITGPI::Acquisition::Seismogram<ValueType>::getTraceType() const
{
    return (seismoType);
}

//! \brief Getter method for reference to coordinates1D vector
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::Seismogram<ValueType>::get1DCoordinates() const
{
    SCAI_ASSERT(coordinates1D.size() == getNumTracesGlobal(), "Size mismatch ");
    return (coordinates1D);
}

//! \brief Getter method for reference to seismogram data
/*!
 *
 * This method returns the DenseMatrix which is used to store the actual seismogram data.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType> &KITGPI::Acquisition::Seismogram<ValueType>::getData()
{
    // SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
    return (data);
}

//! \brief Getter method for const reference to seismogram data
/*!
 *
 * This method returns the DenseMatrix which is used to store the actual seismogram data.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType> const &KITGPI::Acquisition::Seismogram<ValueType>::getData() const
{
    // SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
    return (data);
}

//! \brief Replicate the seismogram data on all processes
/*!
 * Creates a copy of the seismogram data on all processe
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples(new scai::dmemo::NoDistribution(getNumSamples()));
    dmemo::DistributionPtr no_dist_Traces(new scai::dmemo::NoDistribution(getNumTracesGlobal()));

    redistribute(no_dist_Traces, no_dist_numSamples);
}

//! \brief Allocate of the seismogram data
/*!
 * Allocates seismogram based on a given distribution of the traces and the number of samples per trace.
 * The data storage of the seismogram will be distributed according to distTraces. Moreover, the number of
 * local and global traces will be determined based on distTraces.
 *
 \param ctx Context
 \param distTraces Distribution for traces
 \param NT Total number of samples per trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distTraces, IndexType NT)
{

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No Seismogram allocation ");
    SCAI_ASSERT_ERROR(distTraces != NULL, "No valid distribution");

    dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));

    data.setContextPtr(ctx);
    coordinates1D.setContextPtr(ctx);
    inverseAGC.setContextPtr(ctx);

    data.allocate(distTraces, no_dist_NT);
    coordinates1D.allocate(distTraces);
    inverseAGC.allocate(distTraces, no_dist_NT);

    resampleMat = scai::lama::identity<scai::lama::CSRSparseMatrix<ValueType>>(NT); // initialize to identity for no resampling
}

//! \brief Reset of the seismogram data
/*!
 * This method sets the seismogra data to zero. 
 * However, the memory will stay allocated, only the content is overwriten by zeros.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetData()
{
    data.scale(0.0);
}

//! \brief Reset of the seismogram
/*!
 * This method clears the seismogram data and coordinates1D and sets numTraces to zero
 * However, the memory will stay allocated, only the content is overwriten by zeros.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetSeismogram()
{
    data.clear();
    coordinates1D = lama::DenseVector<scai::IndexType>();
    sourceIndex = 0;
    DT = 0.0;
}

//! \brief Redistribute the seismogram data
/*!
 *
 * Redistribution of the seismogram data according to the given distributions.
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::redistribute(scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    if (distSamples == NULL) {
        SCAI_ASSERT_DEBUG(getNumSamples() >= 0, "data has 0 collumns");
        dmemo::DistributionPtr distSamplestmp(new scai::dmemo::NoDistribution(getNumSamples()));
        distSamples = distSamplestmp;
    }

    data.redistribute(distTraces, distSamples);
    coordinates1D.redistribute(distTraces);
}

/*! \brief Getter method for the temporal sampling
 *
 * This method returns the temporal sampling DT in seconds.
 \return DT in seconds
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Seismogram<ValueType>::getDT() const
{
    SCAI_ASSERT_ERROR(DT != 0, "Seismogramm data is not allocated");
    return (DT);
}

/*! \brief Getter method for the number of samples per trace
 *
 *
 \return The number of samples per trace
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumSamples() const
{
    return (data.getNumColumns());
}

/*! \brief Getter method for the number of local traces
 *
 *
 \return The number of local traces on this process
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesLocal() const
{
    return (data.getRowDistribution().getLocalSize());
}

/*! \brief Getter method for the number of global traces
*
*
\return The number of global traces of this seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesGlobal() const
{
    return (data.getNumRows());
}

/*! \brief Getter method for the source Index
*
*
\return The source Index
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getSourceCoordinate() const
{
    return (sourceIndex);
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator*=(ValueType const &rhs)
{
    data *= rhs;

    return *this;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator*=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    data.binaryOp(data, scai::common::BinaryOp::MULT, rhs.data);

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator+(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const
{
    KITGPI::Acquisition::Seismogram<ValueType> result(*this);
    result += rhs;

    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator+=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    data += rhs.data;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator-(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const
{
    KITGPI::Acquisition::Seismogram<ValueType> result(*this);
    result -= rhs;

    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator-=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    data -= rhs.data;

    return *this;
}

/*! \brief Overloading copy assignment operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> &KITGPI::Acquisition::Seismogram<ValueType>::operator=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    //copy rhs with copy constructor to tmp seismogram
    KITGPI::Acquisition::Seismogram<ValueType> tmp(rhs);
    //swap tmp with *this (lhs)
    swap(tmp);

    return *this;
}

template class KITGPI::Acquisition::Seismogram<double>;
template class KITGPI::Acquisition::Seismogram<float>;
