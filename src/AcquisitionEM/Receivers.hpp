#pragma once

#include "AcquisitionGeometry.hpp"
#include "../Acquisition/AcquisitionSettings.hpp"
#include "../Acquisition/suHandler.hpp"
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "Sources.hpp"

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
         * afterwards it will determine the individual #SeismogramTypeEM of each receiver, and based on the #SeismogramTypeEM it will
         * initialize the #Seismogram's of the #SeismogramHandler #receiver.\n\n
         * The #coordinates and #receiver_type vectors contain the coordinates and #SeismogramTypeEM of all receivers, whereas the
         * individual Seismogram of the SeismogramHandler #receiver will hold the receivers seperated based on the #SeismogramTypeEM.
         *
         */
        template <typename ValueType>
        class ReceiversEM : public AcquisitionGeometryEM<ValueType>
        {
          public:
            //! \brief Default constructor
            ReceiversEM(){};

            //! \brief Default destructor
            ~ReceiversEM(){};

            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);

            void init(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber);
            void init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber, scai::IndexType numshots);
            void init(std::vector<receiverSettings> allSettings, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void initWholeSpace(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);

            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings);
            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber);
            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber, scai::IndexType numshots);

          private:
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;
            void acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<receiverSettings> &allSettings, scai::dmemo::DistributionPtr dist_wavefield);
            Acquisition::suHandler<ValueType> su;
        };
    }
}
