#include "Receivers.hpp"

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, config.get<scai::IndexType>("NX"), config.get<scai::IndexType>("NY"), config.get<scai::IndexType>("NZ"), dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<scai::IndexType>("seismoSampleCoeff"));
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(config.get<std::string>("ReceiverFilename") + ".mtx");
    
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, config.get<scai::IndexType>("NX"), config.get<scai::IndexType>("NY"), config.get<scai::IndexType>("NZ"), dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<scai::IndexType>("seismoSampleCoeff"));
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber)
{
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber) + ".mtx");
    
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, config.get<scai::IndexType>("NX"), config.get<scai::IndexType>("NY"), config.get<scai::IndexType>("NZ"), dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<scai::IndexType>("NormalizeTraces"));
    this->getSeismogramHandler().setResampleCoeff(config.get<scai::IndexType>("seismoSampleCoeff"));
}

template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::checkRequiredNumParameter(scai::IndexType numParameterCheck)
{
    /* Check if number of parameters is supported */
    if (numParameterCheck != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }
}

template class KITGPI::Acquisition::Receivers<double>;
template class KITGPI::Acquisition::Receivers<float>;
