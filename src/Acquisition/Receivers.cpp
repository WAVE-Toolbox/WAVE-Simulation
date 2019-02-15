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

    this->setAcquisition(acquisition_temp, modelCoordinates, dist_wavefield, ctx);

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
    scai::lama::DenseMatrix<ValueType> acquisition_temp;

    if (config.get<bool>("initReceiverFromSU")) {
        su.buildAcqMatrixReceiver(config.get<std::string>("ReceiverFilename"), modelCoordinates.getDH());
        acquisition_temp = su.getAcquisition();
    } else
        acquisition_temp.readFromFile(config.get<std::string>("ReceiverFilename") + ".mtx");

    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, modelCoordinates, dist_wavefield, ctx);

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
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    if (config.get<bool>("initReceiverFromSU")) {
        su.buildAcqMatrixReceiver(config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates.getDH());
        acquisition_temp = su.getAcquisition();
    } else
        acquisition_temp.readFromFile(config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber) + ".mtx");
    
    if (config.get<bool>("useReceiversPerShot")) {
        useReceiversPerShot = true;
        shotNr = shotNumber;
    }

    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, modelCoordinates, dist_wavefield, ctx);

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
void KITGPI::Acquisition::Receivers<ValueType>::getAcquisitionMat(Configuration::Configuration const &config, scai::lama::DenseMatrix<ValueType> &acqMat) const
{
    if (config.get<bool>("initSourcesFromSU"))
        acqMat = su.getAcquisition();
    else {
        if (useReceiversPerShot) {
            acqMat.readFromFile(config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNr) + ".mtx"); }
        else {
            acqMat.readFromFile(config.get<std::string>("ReceiverFilename") + ".mtx"); }
    }
}

template class KITGPI::Acquisition::Receivers<double>;
template class KITGPI::Acquisition::Receivers<float>;
