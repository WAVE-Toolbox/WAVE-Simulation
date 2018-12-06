#include "suHandler.hpp"

using namespace scai;

//! \brief Build a Source Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSource(std::string const &filename, ValueType DH)
{
    buildAcqMatrix(filename, DH, &suHandler<ValueType>::buildAcqMatrixSourceComp);
}

//! \brief Build a Receiver Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiver(std::string const &filename, ValueType DH)
{
    buildAcqMatrix(filename, DH, &suHandler<ValueType>::buildAcqMatrixReceiverComp);
}

//! \brief Build a Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 \param buildAcqMat Function pointer to decide if a source or receiver matrix should be build
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrix(std::string const &filename, ValueType DH, buildAcqMatPtr buildAcqMat)
{
    acqMat.clear();
    std::vector<lama::DenseMatrix<ValueType>> acqMatVec(NUM_ELEMENTS_SEISMOGRAMTYPE);

    std::string filenameTmp;

    // read all source files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".SU";
        (this->*buildAcqMat)(filenameTmp, acqMatVec[iComponent], DH);
        nShots[iComponent] = acqMatVec[iComponent].getNumRows();
    }

    // count number of sources to allocate final acquisition matrix
    IndexType numSources = 0;
    for (auto it = acqMatVec.begin(); it != acqMatVec.end(); it++)
        numSources += (*it).getNumRows();

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(numSources, dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(acqMatVec[0].getNumColumns());

    acqMat.allocate(rowDist, colDist);

    // concatenate the acquisition matrices of all components
    auto write_acquisition_HA = hostWriteAccess(acqMat.getLocalStorage().getData());

    IndexType j = 0;
    for (auto it = acqMatVec.begin(); it != acqMatVec.end(); it++) {
        auto read_acquisition_HA = hostReadAccess((*it).getLocalStorage().getData());
        for (IndexType i = 0; i < read_acquisition_HA.size(); i++) {
            write_acquisition_HA[j] = read_acquisition_HA[i];
            j++;
        }
        read_acquisition_HA.release();
    }
    write_acquisition_HA.release();
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSourceComp(std::string const &filename, lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU<ValueType>(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(header.size(), dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(9);

    acqMatTmp.allocate(rowDist, colDist);

    auto write_acquisition_HA = writeAccess(acqMatTmp.getLocalStorage().getData());

    for (int i = 0; i < write_acquisition_HA.size() / 9; i++) {
        thisHeader = header[i];
        write_acquisition_HA[i * 9 + 0] = common::Math::floor<ValueType>(thisHeader.sx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        write_acquisition_HA[i * 9 + 1] = common::Math::floor<ValueType>(thisHeader.sdepth * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH + 0.5);
        write_acquisition_HA[i * 9 + 2] = common::Math::floor<ValueType>(thisHeader.sy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        write_acquisition_HA[i * 9 + 3] = ValueType(component);
        write_acquisition_HA[i * 9 + 4] = 3; // each source signal should be read from file
    }

    write_acquisition_HA.release();
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiverComp(std::string const &filename, lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU<ValueType>(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(header.size(), dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(4);

    acqMatTmp.allocate(rowDist, colDist);

    auto write_acquisition_HA = writeAccess(acqMatTmp.getLocalStorage().getData());

    for (int i = 0; i < write_acquisition_HA.size() / 4; i++) {
        thisHeader = header[i];
        write_acquisition_HA[i * 4 + 0] = common::Math::floor<ValueType>(thisHeader.gx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        write_acquisition_HA[i * 4 + 1] = common::Math::floor<ValueType>(thisHeader.gelev * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH);
        write_acquisition_HA[i * 4 + 2] = common::Math::floor<ValueType>(thisHeader.gy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        write_acquisition_HA[i * 4 + 3] = ValueType(component);
    }

    write_acquisition_HA.release();
}

//! \brief Locate the trace in the SU files based on a shot number
/*!
 \param filename Filename which gets the correct suffix
 \param traceNumber Trace number the shot has in file filename
 \param shotNumber Shot which should be located
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::locateTrace(std::string &filename, scai::IndexType &traceNumber, scai::IndexType shotNumber)
{
    IndexType iComponent = 0;
    do {
        shotNumber -= nShots[iComponent];
        iComponent++;
    } while (shotNumber >= 0);

    filename += "." + std::string(SeismogramTypeString[SeismogramType(iComponent - 1)]) + ".SU";
    traceNumber = nShots[iComponent - 1] + shotNumber;
}

//! \brief Derive the component of the SU file based on its filename
/*!
 \param filename Filename where the component is the suffix
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::suHandler<ValueType>::getComponentFromName(std::string const &filename)
{

    IndexType iTmp = filename.find_last_of('.');
    std::string tmpString = filename.substr(0, iTmp);
    iTmp = tmpString.find_last_of('.');
    tmpString = tmpString.substr(iTmp + 1);

    IndexType component = NUM_ELEMENTS_SEISMOGRAMTYPE + 2;

    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; i++)
        if (tmpString.compare(SeismogramTypeString[SeismogramType(i)]) == 0)
            component = i + 1; // +1 because components in acquisition matrix are defined as p := 1 ...

    SCAI_ASSERT(component < NUM_ELEMENTS_SEISMOGRAMTYPE + 2, "no component found in filename")

    return component;
}

//! \brief Getter member function for the Acquisition Matrix
template <typename ValueType>
lama::DenseMatrix<ValueType> const &KITGPI::Acquisition::suHandler<ValueType>::getAcquisition() const
{
    return (acqMat);
}

//! \brief Getter member function for the Acquisition Matrix
/*!
 \param acqRowMat Acquisition Matrix
 \param shotNumber Shot number to get
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::getAcquisitionRow(lama::DenseMatrix<ValueType> &acqRowMat, IndexType shotNumber) const
{
    scai::lama::DenseVector<ValueType> acqRow;
    acqMat.getRow(acqRow, shotNumber);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(1, dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(9);

    acqRowMat.allocate(rowDist, colDist);
    acqRowMat.setRow(acqRow, 0, common::BinaryOp::COPY);
}

template class KITGPI::Acquisition::suHandler<double>;
template class KITGPI::Acquisition::suHandler<float>;