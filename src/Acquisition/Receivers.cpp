#include "Receivers.hpp"

/*! \brief Constructor based on the configuration class and the distribution of the wavefields. This constructor will read the acquistion from the Receiverfile
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
KITGPI::Acquisition::Receivers<ValueType>::Receivers(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{

    /* Read acquisition matrix */
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(config.get<std::string>("ReceiverFilename"));
    init(acquisition_temp, config, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(acquisition_temp, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
}

template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::checkRequiredNumParameter(IndexType numParameterCheck)
{
    /* Check if number of parameters is supported */
    if (numParameterCheck != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }
}

template class KITGPI::Acquisition::Receivers<double>;
template class KITGPI::Acquisition::Receivers<float>;
