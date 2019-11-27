#include "suHandler.hpp"
#include "../Common/HostPrint.hpp"
#include "../IO/SUIO.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CollectiveFile.hpp>
using namespace scai;

//! \brief Build a Source Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::readAllSettingsFromSU(std::vector<sourceSettings<ValueType>> &allSettings, std::string const &filename, ValueType DH)
{

    allSettings.clear();
    std::vector<sourceSettings<ValueType>> sourceSettingsVecTmp;

    std::string filenameTmp;

    // read all source files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".su";
        readAllSettingsFromSUComp(filenameTmp, sourceSettingsVecTmp, DH);
        allSettings.insert(allSettings.end(), sourceSettingsVecTmp.begin(), sourceSettingsVecTmp.end());
    }
}

//! \brief Build a Receiver Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::readAllSettingsFromSU(std::vector<receiverSettings> &allSettings, std::string const &filename, ValueType DH)
{
    allSettings.clear();
    std::vector<receiverSettings> receiverSettingsVecTmp;

    std::string filenameTmp;
    IndexType counter = 0;

    // read all receiver files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".su";
        try {
            readAllSettingsFromSUComp(filenameTmp, receiverSettingsVecTmp, DH);
        } catch (scai::common::Exception &e) {
            counter++;
            if (counter == 4) {
                COMMON_THROWEXCEPTION("No file with name: " << filename << ".'comp'.su could be read");
            }
        }
        allSettings.insert(allSettings.end(), receiverSettingsVecTmp.begin(), receiverSettingsVecTmp.end());
    }
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param sourceSettingsVec Vector of source settings
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::readAllSettingsFromSUComp(std::string const &filename, std::vector<sourceSettings<ValueType>> &sourceSettingsVec, ValueType DH)
{
    sourceSettingsVec.clear();
    std::vector<Segy> header;
    SUIO::readHeaderSU<ValueType>(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);
    //sourceSettingsVec.reserve(header.size());

    //sourceSettings<ValueType> sourceSettingsTmp;
    sourceSettingsVec.resize(header.size());
    for (unsigned int i = 0; i < header.size(); i++) {
        thisHeader = header[i];
        sourceSettingsVec[i].sourceCoords.x = static_cast<IndexType>(thisHeader.sx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        sourceSettingsVec[i].sourceCoords.y = static_cast<IndexType>(thisHeader.sdepth * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH + 0.5);
        sourceSettingsVec[i].sourceCoords.z = static_cast<IndexType>(thisHeader.sy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        sourceSettingsVec[i].sourceType = component;
        sourceSettingsVec[i].waveletType = 3; // each source signal should be read from file
        //    sourceSettingsVec.push_back(sourceSettingsTmp);
    }
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::readAllSettingsFromSUComp(std::string const &filename, std::vector<receiverSettings> &receiverSettingsVec, ValueType DH)
{
    receiverSettingsVec.clear();
    std::vector<Segy> header;
    SUIO::readHeaderSU<ValueType>(filename, header);
    Segy thisHeader;
    IndexType component = getComponentFromName(filename);

    receiverSettingsVec.resize(header.size());

    for (unsigned int i = 0; i < header.size(); i++) {
        thisHeader = header[i];
        receiverSettingsVec[i].receiverCoords.x = static_cast<IndexType>(thisHeader.gx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        receiverSettingsVec[i].receiverCoords.y = static_cast<IndexType>(thisHeader.gelev * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH);
        receiverSettingsVec[i].receiverCoords.z = static_cast<IndexType>(thisHeader.gy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        receiverSettingsVec[i].receiverType = component;
    }
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

template class KITGPI::Acquisition::suHandler<double>;
template class KITGPI::Acquisition::suHandler<float>;
