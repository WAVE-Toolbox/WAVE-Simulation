#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>
#include <vector>

#include "AcquisitionGeometry.hpp"
#include "../Acquisition/AcquisitionSettings.hpp"
#include "../Acquisition/SourceSignal/all.hpp"
#include "../Acquisition/suHandler.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        template <typename ValueType>
        class AcquisitionGeometryEM;

        template <typename ValueType>
        class SeismogramEM;

        template <typename ValueType>
        class SeismogramHandlerEM;

        //! \brief Handling of sources
        /*!
         * This class accounts for the handling of seismic sources.
         * It provides the reading from the source acquisition from file, the distribution of the sources and the generation of synthetic signals.
         */
        template <typename ValueType>
        class SourcesEM : public AcquisitionGeometryEM<ValueType>
        {

          public:
            //! \brief Default destructor
            ~SourcesEM(){};
            void init(std::vector<sourceSettings<ValueType>> allSettings, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseMatrix<ValueType> &signalMatrix);

            void generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, std::vector<scai::IndexType>);

            void getAcquisitionSettings(Configuration::Configuration const &config, std::vector<sourceSettings<ValueType>> &allSettings);
            scai::lama::DenseMatrix<ValueType> getsourcesignal();

          private:
            SeismogramEM<ValueType> signals; //!< Source signals
            Acquisition::suHandler<ValueType> su;

            /* Requiered acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_type; //!< Type of wavelet: 1==Synthetic

            /* Optional acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_shape; //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            scai::lama::DenseVector<ValueType> wavelet_fc;          //!< Center frequency of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_amp;         //!< Amplitude of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_tshift;      //!< Time shift of synthetic wavelet

            bool wavelet_type_flag_2 = false; // flag if wavelet type 2 is used
            bool wavelet_type_flag_3 = false; // flag if wavelet type 3 is used

            void initOptionalAcquisitionParameter(scai::IndexType numParameter, scai::IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx) override;
            void initOptionalAcquisitionParameter(std::vector<sourceSettings<ValueType>> allSettings, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx);
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;
            void acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<sourceSettings<ValueType>> &allSettings, scai::dmemo::DistributionPtr dist_wavefield);

            void copySignalsToSeismogramHandler();

            void allocateSeismogram(scai::IndexType NT, scai::dmemo::DistributionPtr dist_traces, scai::hmemo::ContextPtr ctx);
            void generateSyntheticSignal(scai::IndexType SourceLocal, scai::IndexType NT, ValueType DT);
            void readSignalFromFile(Configuration::Configuration const &config, scai::IndexType SourceLocal, scai::IndexType rowNumber);
        };
    }
}
