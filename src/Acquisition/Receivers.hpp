#pragma once

#include "AcquisitionGeometry.hpp"
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include "suHandler.hpp"

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

            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            
            void init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber);
            
            void getAcquisitionMat(Configuration::Configuration const &config, scai::lama::DenseMatrix<ValueType> &acqMat) const;

          private:
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;
            suHandler<ValueType> su;
        };
    }
}
