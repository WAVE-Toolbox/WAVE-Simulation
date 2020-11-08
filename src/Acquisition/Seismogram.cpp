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
    data.swap(rhs.data);
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

        std::string filenameTmp = filename + "." + SeismogramTypeString[getTraceType()];

        switch (seismogramFormat) {
        case 4:
            SUIO::writeSU(filenameTmp, dataResample, coordinates1D, outputDT, sourceIndex, modelCoordinates);
            break;
        default:
            IO::writeMatrix(dataResample, filenameTmp, seismogramFormat);
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
    switch (seismogramFormat) {
    case 4:
        SUIO::readDataSU(filenameTmp, data, this->getNumSamples(), this->getNumTracesGlobal());
        break;
    default:
        IO::readMatrix(data, filenameTmp, seismogramFormat);
        break;
    }
}

//! \brief Normalize the seismogram-traces
/*!
 *
 * This methode normalized the traces of the seismogram after the time stepping.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::normalizeTrace()
{

    scai::hmemo::HArray<ValueType> tempRow;
    for (IndexType i = 0; i < getNumTracesLocal(); i++) {
        data.getLocalStorage().getRow(tempRow, i);
        ValueType tempMax = scai::utilskernel::HArrayUtils::max(tempRow);
        if (tempMax!=0) {
            scai::utilskernel::HArrayUtils::setScalar(tempRow, tempMax, scai::common::BinaryOp::DIVIDE);
            data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
        }
    }
}

//! \brief Set killed traces to zero
/*!
 *
 * This methode sets the killed traces from the preprocessing to zero
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::killTrace()
{
    if (data.getNumRows() != 0) {
        std::vector<scai::IndexType> deadTraces = {0, 9, 11};
        scai::lama::DenseVector<ValueType> tempRow(data.getNumColumns(), 0);
        for (IndexType i = 0; i < 3; i++) {
            IndexType deadTraceInd = deadTraces[i];
            data.setRow(tempRow, deadTraceInd, scai::common::BinaryOp::COPY);
        }
    }

}

//! \brief Integrate the seismogram-traces
/*!
 *
 * This methode integrate the traces of the seismogram.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::integrateTraces()
{

    scai::hmemo::HArray<ValueType> tempRow;
    for (IndexType i = 0; i < getNumTracesLocal(); i++) {
        data.getLocalStorage().getRow(tempRow, i);
        for (IndexType j = 0; j < tempRow.size() - 1; j++) {
            tempRow[j + 1] = tempRow[j + 1] * DT + tempRow[j];
        }
        data.getLocalStorage().setRow(tempRow, i, scai::common::BinaryOp::COPY);
    }
}

//! \brief Filter the seismogram-traces
/*!
 *
 * This methode filters the traces of the seismogram.
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

//! \brief Setter methode to set matrix for resampling this seismogram.
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

    data.allocate(distTraces, no_dist_NT);
    coordinates1D.allocate(distTraces);

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
