#pragma once

#include <scai/lama.hpp>

#include "AcquisitionGeometry.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        /*! \brief Handling of receivers
         *
         * This class accounts for the handling of seismic receivers.
         * It provides the reading of the receivers acquisition from file, the distribution of the receivers and the collection of the seismograms.
         *
         *
         * This class will first read-in the global receiver configuration from a receiver configuration file and
         * afterwards it will determine the individual #SeismogramType of each receiver, and based on the #SeismogramType it will
         * initialize the #Seismogram's of the #SeismogramHandler #receiver.\n\n
         * The #coordinates and #receiver_type vectors contain the coordinates and #SeismogramType of all receivers, whereas the
         * individual Seismogram of the SeismogramHandler #receiver will hold the receivers seperated based on the #SeismogramType.
         *
         */
        template <typename ValueType>
        class Receivers : public AcquisitionGeometry<ValueType>
        {
          public:
            //! \brief Default constructor
            Receivers(){};

            explicit Receivers(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);

            //! \brief Default destructor
            ~Receivers(){};

            void init(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield);

          private:
            void checkRequiredNumParameter(IndexType numParameterCheck);
        };
    }
}

/*! \brief Constructor based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
KITGPI::Acquisition::Receivers<ValueType>::Receivers(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
{
    init(config, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist_wavefield)
{
    IndexType getNT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->readAcquisitionFromFile(config.get<std::string>("ReceiverFilename"), config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), dist_wavefield, ctx);

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
