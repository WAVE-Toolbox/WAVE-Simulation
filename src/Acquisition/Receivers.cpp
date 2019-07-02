#include "Receivers.hpp"

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    std::vector<receiverSettings> allSettings;
    acqMat2settings(acquisition_temp, allSettings, dist_wavefield);
    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<ValueType>("seismoDT") / config.get<ValueType>("DT"));
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    std::vector<receiverSettings> allSettings;

    if (config.get<bool>("initReceiverFromSU")) {
        su.buildAcqMatrixReceiver(config.get<std::string>("ReceiverFilename"), modelCoordinates.getDH());
        allSettings = su.getReceiverSettingsVec();
    } else
        readAllSettings(allSettings, std::string(config.get<std::string>("ReceiverFilename") + ".txt"));

    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<ValueType>("seismoDT") / config.get<ValueType>("DT"));
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber)
{
    std::vector<receiverSettings> allSettings;
    if (config.get<bool>("initReceiverFromSU")) {
        su.buildAcqMatrixReceiver(config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates.getDH());
        allSettings = su.getReceiverSettingsVec();
    } else
        readAllSettings(allSettings, config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber) + ".txt");

    if (config.get<bool>("useReceiversPerShot")) {
        useReceiversPerShot = true;
        shotNr = shotNumber;
    }

    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<ValueType>("seismoDT") / config.get<ValueType>("DT"));
}

template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::checkRequiredNumParameter(scai::IndexType numParameterCheck)
{
    /* Check if number of parameters is supported */
    if (numParameterCheck != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by MTX or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getAcquisitionMat(Configuration::Configuration const &config, std::vector<receiverSettings> &allReceiverSettings)
{
    if (config.get<bool>("initSourcesFromSU"))
        allReceiverSettings = su.getReceiverSettingsVec();
    else {
        if (useReceiversPerShot) {
            readAllSettings(allReceiverSettings, config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNr) + ".txt");
        } else {
            readAllSettings(allReceiverSettings, config.get<std::string>("ReceiverFilename") + ".txt");
        }
    }
}

template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<receiverSettings> &allSettings, scai::dmemo::DistributionPtr dist_wavefield)
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

template class KITGPI::Acquisition::Receivers<double>;
template class KITGPI::Acquisition::Receivers<float>;
