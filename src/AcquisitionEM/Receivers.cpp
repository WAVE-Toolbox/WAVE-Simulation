#include "Receivers.hpp"
#include "../CheckParameter/CheckParameter.hpp"
#include "../IO/IO.hpp"

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param allSettings allSettings for receivers
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::init(std::vector<receiverSettings> allSettings, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{    
    /*reset seismograms. This is necessary when init will be called multiple times*/
//     this->getSeismogramHandler().resetSeismograms();
    
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setSeismoDT(config.get<ValueType>("seismoDT"));
    this->getSeismogramHandler().setEnvelopTrace(config.get<ValueType>("envelopTraces"));
}

template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    std::vector<receiverSettings> allSettings;
    acqMat2settings(acquisition_temp, allSettings, dist_wavefield);
    
    CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::init(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    std::vector<receiverSettings> allSettings;
    
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    if (useStreamConfig) {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        std::vector<receiverSettings> allSettingsBig;
        getAcquisitionSettings(configBig, allSettingsBig);
        
        scai::IndexType shotIndTrue = 0;
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;  
        Acquisition::Sources<ValueType> sources;
        sources.getAcquisitionSettings(configBig, sourceSettingsBig);
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
        Acquisition::getSettingsPerShot<ValueType>(allSettings, allSettingsBig, cutCoordinates, shotIndTrue);
    } else {
        getAcquisitionSettings(config, allSettings);
    }

    CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber)
{
    std::vector<receiverSettings> allSettings;
    
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    if (useStreamConfig) {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        std::vector<receiverSettings> allSettingsBig;
        getAcquisitionSettings(configBig, allSettingsBig, shotNumber);
        
        scai::IndexType shotIndTrue = 0;
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;  
        Acquisition::Sources<ValueType> sources;
        sources.getAcquisitionSettings(configBig, sourceSettingsBig);
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
        Acquisition::getSettingsPerShot<ValueType>(allSettings, allSettingsBig, cutCoordinates, shotIndTrue);
    } else {
        getAcquisitionSettings(config, allSettings, shotNumber);
    }
    
    try {
        CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
    } catch (scai::common::Exception &e) {
        HOST_PRINT(dist_wavefield->getCommunicatorPtr(), "Message from receiver.cpp init: Error while reading receiver settings for shot " << shotNumber);
        COMMON_THROWEXCEPTION("e.what()");
    }
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber, scai::IndexType numshots)
{
    std::vector<receiverSettings> allSettings;
    
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    if (useStreamConfig) {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        std::vector<receiverSettings> allSettingsBig;
        getAcquisitionSettings(configBig, allSettingsBig, shotNumber, numshots);
        
        scai::IndexType shotIndTrue = 0;
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;  
        Acquisition::Sources<ValueType> sources;
        sources.getAcquisitionSettings(configBig, sourceSettingsBig);
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
        Acquisition::getuniqueShotInd(shotIndTrue, uniqueShotNos, shotNumber);
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        Acquisition::getCutCoord(cutCoordinates, sourceSettingsBig);
        Acquisition::getSettingsPerShot<ValueType>(allSettings, allSettingsBig, cutCoordinates, shotIndTrue);
    } else {
        getAcquisitionSettings(config, allSettings, shotNumber, numshots);
    }
    
    try {
        CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
    } catch (scai::common::Exception &e) {
        HOST_PRINT(dist_wavefield->getCommunicatorPtr(), "Message from receiver.cpp init: Error while reading receiver settings for shot " << shotNumber);
        COMMON_THROWEXCEPTION("e.what()");
    }
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::checkRequiredNumParameter(scai::IndexType numParameterCheck)
{
    /* Check if number of parameters is supported */
    if (numParameterCheck != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings)
{
    if (config.get<bool>("initReceiverFromSU")) {
        su.readAllSettingsFromSU(allSettings, config.get<std::string>("ReceiverFilename"), config.get<ValueType>("DH"));
    } else if (config.get<scai::IndexType>("useReceiversPerShot") == 0) {
        readAllSettings(allSettings, config.get<std::string>("ReceiverFilename") + ".txt");
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber)
{
    if (config.get<scai::IndexType>("useReceiversPerShot") == 1) {
        if (config.get<bool>("initReceiverFromSU")) {
            su.readAllSettingsFromSU(allSettings, config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber), config.get<ValueType>("DH"));
        } else {
            readAllSettings(allSettings, config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber) + ".txt");
        }
    } 
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber, scai::IndexType numshots)
{
    if (config.get<scai::IndexType>("useReceiversPerShot") == 2) {
        std::vector<receiverSettings> allSettingsTemp;
        allSettings.clear();
        if (config.get<bool>("initReceiverFromSU")) {
            su.readAllSettingsFromSU(allSettingsTemp, config.get<std::string>("ReceiverFilename"), config.get<ValueType>("DH"));
        } else {
            readAllSettings(allSettingsTemp, std::string(config.get<std::string>("ReceiverFilename") + ".txt"));
        }
        std::string filenameTmp = config.get<std::string>("ReceiverFilename") + ".mark";
        scai::lama::CSRSparseMatrix<ValueType> receiverMarkMatrix;
        receiverMarkMatrix.allocate(numshots, allSettingsTemp.size()+1);
        IO::readMatrix<ValueType>(receiverMarkMatrix, filenameTmp, 1); // .mtx file
        scai::lama::SparseVector<ValueType> receiverMarkVector;
        for (scai::IndexType shotInd = 0; shotInd < numshots; shotInd++) {
            receiverMarkMatrix.getRow(receiverMarkVector, shotInd);
            if (receiverMarkVector[0] == shotNumber) {
                break;
            } else if (shotInd == numshots - 1) {
                SCAI_ASSERT(receiverMarkVector[0] == shotNumber, "receiverMarkVector[0] != shotNumber");
            }
        }
                    
        for (scai::IndexType iMark = 0; iMark < receiverMarkVector.size()-1; iMark++) {
            if (receiverMarkVector[iMark+1] != 0) {
                allSettings.push_back(allSettingsTemp[iMark]);
            }
        }
    }
}

template <typename ValueType>
void KITGPI::Acquisition::ReceiversEM<ValueType>::acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<receiverSettings> &allSettings, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::IndexType nRows = acqMat.getNumRows();
    allSettings.reserve(nRows);
    if (dist_wavefield->getCommunicator().getRank() == 0) {
        auto read_acquisition_temp_HA = hostReadAccess(acqMat.getLocalStorage().getData());
        for (scai::IndexType row = 0; row < nRows; row++) {
            allSettings[row].receiverCoords.x = read_acquisition_temp_HA[row * 9 + 0];
            allSettings[row].receiverCoords.y = read_acquisition_temp_HA[row * 9 + 1];
            allSettings[row].receiverCoords.z = read_acquisition_temp_HA[row * 9 + 2];
            allSettings[row].receiverType = read_acquisition_temp_HA[row * 9 + 3];
        }
    }
}

template class KITGPI::Acquisition::ReceiversEM<double>;
template class KITGPI::Acquisition::ReceiversEM<float>;
