#pragma once

#include "AcquisitionGeometry.hpp"
#include "AcquisitionSettings.hpp"
#include "suHandler.hpp"
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

            //! \brief Default destructor
            ~Receivers(){};

            void init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);

            void init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber, std::vector<sourceSettings<ValueType>> sourceSettingsEncode);
            void init(std::vector<receiverSettings> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void initWholeSpace(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseVector<scai::IndexType> const &seismogramTypes);

            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings);
            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber);
            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber, scai::IndexType numshots, std::vector<IndexType> shotIndIncr, std::vector<sourceSettings<ValueType>> sourceSettingsEncode);
            void writeReceiverMark(std::string filename);
            void decode(Configuration::Configuration const &config, std::string const &filename, scai::IndexType shotNumber, std::vector<sourceSettings<ValueType>> sourceSettingsEncode);

            using AcquisitionGeometry<ValueType>::isSeismic;
            
          private:
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;
            void acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<receiverSettings> &allSettings, scai::dmemo::DistributionPtr dist_wavefield);
            suHandler<ValueType> su;
            scai::lama::SparseVector<ValueType> receiverMarkVector;
        };
    }
}
