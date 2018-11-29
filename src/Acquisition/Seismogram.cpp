#include "Seismogram.hpp"
using namespace scai;

//! \brief copy constructor
/*!
 *
 \param rhs seismogram to copy
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>::Seismogram(const Seismogram &rhs)
{
    //   std::cout<< "copy seismogram" << std::endl;
    numSamples = rhs.numSamples;
    numTracesGlobal = rhs.numTracesGlobal;
    numTracesLocal = rhs.numTracesLocal;
    DT = rhs.DT;
    type = rhs.type;
    coordinates = rhs.coordinates;
    sourceCoordinate = rhs.sourceCoordinate;
    data = rhs.data;
    resampleMatRight = rhs.resampleMatRight;
    resampleMatLeft = rhs.resampleMatLeft;
    resampleVec = rhs.resampleVec;
    //   std::cout<< "copy seismogram data " << data << std::endl;
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
    std::swap(numSamples, rhs.numSamples);
    std::swap(numTracesGlobal, rhs.numTracesGlobal);
    std::swap(numTracesLocal, rhs.numTracesLocal);
    std::swap(DT, rhs.DT);
    std::swap(type, rhs.type);
    std::swap(coordinates, rhs.coordinates);
    std::swap(sourceCoordinate, rhs.sourceCoordinate);
    data.swap(rhs.data);
}

//! \brief Adding ending to the seismogram-filename-string
/*!
 *
 * This member function adds the #SeismogramType to the filname ending.
 \param filename Filename of output
 */
template <typename ValueType>
std::string KITGPI::Acquisition::Seismogram<ValueType>::addSeismogramTypeToName(std::string const &filename) const
{
    SCAI_ASSERT_DEBUG(type >= 0 && type <= (NUM_ELEMENTS_SEISMOGRAMTYPE - 1), "Wrong Trace Type: " << getTraceType());
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::P], "p") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VX], "vx") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VY], "vy") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VZ], "vz") == 0, "Error in mapping of SeismogramType to std::string");

    std::size_t found = filename.find_last_of(".");
    std::string beforeEnding = filename.substr(0, found);
    std::string afterEnding = filename.substr(found);
    std::string traceTypeString = SeismogramTypeString[getTraceType()];
    return (beforeEnding + "." + traceTypeString + afterEnding);
}

//! \brief Setter method for the context ptr
/*!
 *
 * This method sets the Context to the coordinates and to the seismogram data.
 \param ctx Set ContextPtr
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setContextPtr(scai::hmemo::ContextPtr ctx)
{
    data.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);
}

//! \brief  Read the seismogram from disk
/*!
 *
 * This method writes the seismogram data to disk. It uses the Configuration class to determine the filename and some requiered header information.\n
 * The output format is determined by the input parameter `SeismogramFormat`, which will be requested to the Configuration class. \n
 * **Supported Formats:**\n
 * 1. MTX: MatrixMaker format
 * 2. SU: SeismicUnix format
 \param config Configuration class which is used to determine the filename and header information
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::write(Configuration::Configuration const &config, std::string const &filename) const
{
    switch (config.get<IndexType>("SeismogramFormat")) {
    case 1:
        writeToFileRaw(filename + ".mtx");
        break;
    case 2:
        writeToFileSU(filename + ".SU", config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        break;
    default:
        COMMON_THROWEXCEPTION(" Unkown SeismogramFormat ")
        break;
    }
}

//! \brief Read the seismogram from disk
/*!
 *
 * This method reads the seismogram data from disk. It uses the Configuration class to determine the format.\n
 * **Supported Formats:**\n
 * 1. MTX: MatrixMaker format
 * 2. SU: SeismicUnix format
 \param config Configuration class which is used to determine the filename and header information
 \param filename Filename to read seismogram
 \param copyDist Boolean: 0 = read data undistributed (default), data is replicated on each process // 1 = read data with existing distribution of data
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::read(Configuration::Configuration const &config, std::string const &filename, bool copyDist)
{
    switch (config.get<IndexType>("SeismogramFormat")) {
    case 1:
        readFromFileRaw(filename + ".mtx", copyDist);
        break;
    case 2:
        readFromFileSU(filename + ".SU", copyDist);
        break;
    default:
        COMMON_THROWEXCEPTION(" Unkown SeismogramFormat ")
        break;
    }
}

//! \brief Read the seismogram from disk
/*!
 *
 * This method reads the seismogram data from disk. It uses the Configuration class to determine the format.\n
 * **Supported Formats:**\n
 * 1. MTX: MatrixMaker format
 * 2. SU: SeismicUnix format
 \param config Configuration class which is used to determine the filename and header information
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::read(Configuration::Configuration const &config, std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    switch (config.get<IndexType>("SeismogramFormat")) {
    case 1:
        readFromFileRaw(filename + ".mtx", distTraces, distSamples);
        break;
    case 2:
        readFromFileSU(filename + ".SU", distTraces, distSamples);
        break;
    default:
        COMMON_THROWEXCEPTION(" Unkown SeismogramFormat ")
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
    if (normalizeTraces == 1) {

        SCAI_ASSERT(data.getNumRows() == numTracesGlobal, " Size of matrix is not matching with number of traces. ");

        scai::lama::DenseVector<ValueType> tempRow;
        ValueType tempMax;
        ValueType tempInverseMax;

        for (IndexType i = 0; i < numTracesGlobal; i++) {
            tempMax = 0.0;
            data.getRow(tempRow, i);
            tempMax = tempRow.max();
            tempInverseMax = 1 / tempMax;
            tempRow *= tempInverseMax;
            data.setRow(tempRow, i, scai::common::BinaryOp::COPY);
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
    std::cout << "integrate traces, data = " << data << std::endl;

    SCAI_ASSERT(data.getNumRows() == numTracesGlobal, " Size of matrix is not matching with number of traces. ");

    scai::lama::DenseVector<ValueType> tempRow;

    for (IndexType i = 0; i < numTracesGlobal; i++) {
        data.getRow(tempRow, i);
        for (IndexType j = 0; j < tempRow.size() - 1; j++) {
            tempRow[j + 1] = tempRow[j + 1] * DT + tempRow[j];
        }
        data.setRow(tempRow, i, scai::common::BinaryOp::COPY);
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

//! \brief Setter method for the source coordinate
/*!
 *
 * This method sets the source coordinate to this class. The source coordinate will be used for calculation of the header information (e.g. Offset), mainly during seismogram output to disk.
 \param sourceCoord Source coordinate in 1-D format
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    SCAI_ASSERT_DEBUG(sourceCoord >= 0, "sourceCoord is not valid");
    sourceCoordinate = sourceCoord;
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
    type = trace;
};

//! \brief Setter function for the coordinates of the traces
/*!
 *
 * This method will set the coordinates of the traces to this class. 
 * The size of the coordinate vector has to be equal to the number of global traces.
 \param coord DenseVector with coordinates
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setCoordinates(scai::lama::DenseVector<IndexType> const &coord)
{
    SCAI_ASSERT_ERROR(coord.size() == numTracesGlobal, "Given traceType vector has wrong format");
    coordinates = coord;
};

//! \brief Setter methode to set Index for trace-normalization.
/*!
 *
 * This method sets the index for trace-normalization.
 \param normalizeTrace Index for trace-normalization which will normalize the seismogram traces
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setNormalizeTraces(IndexType normalize)
{
    SCAI_ASSERT(normalize >= 0 && normalize <= 1, " Index has to be 1 or 0 ");
    normalizeTraces = normalize;
}

//! \brief Setter methode to set matrix for resampling this seismogram.
/*!
 \param rMat Resampling matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setResampleCoeff(ValueType resampleCoeff)
{
    if (this->getNumSamples() != 0) {
        Common::calcResampleMat<ValueType>(resampleMatLeft, numSamples, resampleCoeff, 0);
        if (scai::common::Math::floor<ValueType>(resampleCoeff) != resampleCoeff) {
            Common::calcResampleMat<ValueType>(resampleMatRight, numSamples, resampleCoeff, 1);
            Common::calcResampleVec<ValueType>(resampleVec, numSamples, resampleCoeff);
        }
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
    return (type);
}

//! \brief Getter method for reference to coordinates vector
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::Seismogram<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_DEBUG(coordinates.size() == numTracesGlobal, "Size mismatch ");
    return (coordinates);
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
    SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
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
    SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
    return (data);
}

//! \brief Replicate the seismogram data on all processes
/*!
 * Creates a copy of the seismogram data on all processe
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples(new scai::dmemo::NoDistribution(numSamples));
    dmemo::DistributionPtr no_dist_Traces(new scai::dmemo::NoDistribution(numTracesGlobal));

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
//       std::cout << "Seismogram allocate dist = " << *distTraces << " x NT = " << NT << std::endl;

    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No Seismogram allocation ");
    SCAI_ASSERT_ERROR(distTraces != NULL, "No valid distribution");

    numSamples = NT;
    numTracesGlobal = distTraces->getGlobalSize();
    numTracesLocal = distTraces->getLocalSize();

    dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));

    data.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);

    data.allocate(distTraces, no_dist_NT);
    coordinates.allocate(distTraces);
    
    resampleMatLeft = scai::lama::identity<scai::lama::CSRSparseMatrix<ValueType>>(NT); // initialize to identity for no resampling
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
 * This method clears the seismogram data and coordinates and sets numTraces to zero
 * However, the memory will stay allocated, only the content is overwriten by zeros.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetSeismogram()
{
    numTracesGlobal = 0;
    numTracesLocal = 0;
    data.clear();
    coordinates = lama::DenseVector<scai::IndexType>();
    sourceCoordinate = 0;
    numSamples = 0;
    DT = 0.0;
    normalizeTraces = 0;
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
        SCAI_ASSERT_DEBUG(numSamples >= 0, "numSamples not set");
        dmemo::DistributionPtr distSamplestmp(new scai::dmemo::NoDistribution(numSamples));
        distSamples = distSamplestmp;
    }

    data.redistribute(distTraces, distSamples);
    coordinates.redistribute(distTraces);
}

//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param copyDist Boolean: 0 = read data undistributed (default), data is replicated on each process // 1 = read data with existing distribution of data
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileRaw(std::string const &filename, bool copyDist)
{
    scai::dmemo::DistributionPtr distTraces;
    scai::dmemo::DistributionPtr distSamples;

    if (copyDist == 1) {
        distTraces = data.getRowDistributionPtr();
        distSamples = data.getColDistributionPtr();
    }

    data.readFromFile(addSeismogramTypeToName(filename));
    IndexType nrow_temp = data.getNumRows();
    IndexType ncolumn_temp = data.getNumColumns();

    numSamples = ncolumn_temp;
    numTracesGlobal = nrow_temp;

    if (copyDist == 0) {
        replicate();
    } else if (copyDist == 1) {
        redistribute(distTraces, distSamples);
    }
}

//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param copyDist Boolean: 0 = read data undistributed (default), data is replicated on each process // 1 = read data with existing distribution of data
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileSU(std::string const &filename, bool copyDist)
{   
    scai::dmemo::DistributionPtr distTraces;
    scai::dmemo::DistributionPtr distSamples;

    if (copyDist == 1) {
        distTraces = data.getRowDistributionPtr();
        distSamples = data.getColDistributionPtr();
    }

    std::string tempstr = addSeismogramTypeToName(filename);
    su.readDataSU(tempstr, data, this->getNumSamples(), this->getNumTracesGlobal());
    
    if (copyDist == 0) {
        replicate();
    } else if (copyDist == 1) {
        redistribute(distTraces, distSamples);
    }
}

//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    data.readFromFile(addSeismogramTypeToName(filename));
    IndexType nrow_temp = data.getNumRows();
    IndexType ncolumn_temp = data.getNumColumns();

    numSamples = ncolumn_temp;
    numTracesGlobal = nrow_temp;

    if (distTraces == NULL && distSamples == NULL) {
        replicate();
    } else {
        redistribute(distTraces, distSamples);
    }
}

//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileSU(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{  
    std::string tempstr = addSeismogramTypeToName(filename);
    su.readDataSU(tempstr, data, this->getNumSamples(), this->getNumTracesGlobal());
      
    if (distTraces == NULL && distSamples == NULL) {
        replicate();
    } else {
        redistribute(distTraces, distSamples);
    }
}

//! \brief Write a seismogram to disk without header
/*!
 *
 \param filename Filename to write seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileRaw(std::string const &filename) const
{
    if (data.getNumValues() > 0) {
        scai::lama::DenseMatrix<ValueType> dataResample;
        dataResample = data * resampleMatLeft;
        
        if (resampleMatRight.getNumColumns() > 0) {
            scai::lama::DenseMatrix<ValueType> dataResampleInt;
            dataResampleInt = data * resampleMatRight;
            dataResampleInt -= dataResample;
            dataResampleInt.scaleColumns(resampleVec);
            dataResample += dataResampleInt;
        }
        dataResample.writeToFile(addSeismogramTypeToName(filename));
    }
}

/*! \brief Getter method for reference normalization index
 *
 *
 \return NormalizeTraces Index 
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNormalizeTraces() const
{
    return (normalizeTraces);
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
    return (numSamples);
}

/*! \brief Getter method for the number of local traces
 *
 *
 \return The number of local traces on this process
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesLocal() const
{
    return (numTracesLocal);
}

/*! \brief Getter method for the number of global traces
*
*
\return The number of global traces of this seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesGlobal() const
{
    return (numTracesGlobal);
}

/*! \brief Getter method for the 1D source coordinate
*
*
\return The 1D source coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getSourceCoordinate() const
{
    return (sourceCoordinate);
}

//! \brief Write a seismogram to disk in Seismic Unix (SEG-Y) format
/*!
 *
 * This method writes the seismogram in the Seismic Unix format to disk.
 * Some header information will be calculated based on the input parameters and will be included in the seismic unix file.
 \param filename Filename to write seismogram in Seismic Unix (SEG-Y) format
 \param NX Number of grid points in X direction
 \param NY Number of grid points in Y direction
 \param NZ Number of grid points in Z direction
 \param DH Length of space step in meter
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileSU(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, ValueType DH) const
{
    if (data.getNumValues() > 0) {
        std::string tempstr = addSeismogramTypeToName(filename);
        su.writeSU(tempstr, data, coordinates, DT, sourceCoordinate, NX, NY, NZ, DH);
    }
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
