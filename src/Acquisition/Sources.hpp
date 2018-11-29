#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "AcquisitionGeometry.hpp"
#include "SourceSignal/all.hpp"
#include "suHandler.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        template <typename ValueType>
        class AcquisitionGeometry;

        template <typename ValueType>
        class Seismogram;

        template <typename ValueType>
        class SeismogramHandler;

        //! \brief Handling of sources
        /*!
         * This class accounts for the handling of seismic sources.
         * It provides the reading from the source acquisition from file, the distribution of the sources and the generation of synthetic signals.
         */
        template <typename ValueType>
        class Sources : public AcquisitionGeometry<ValueType>
        {

          public:
            //! \brief Default constructor
            Sources(){};

            explicit Sources(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            explicit Sources(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber);

            //! \brief Default destructor
            ~Sources(){};
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber);

            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber);
            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseMatrix<ValueType> &signalMatrix);

            void generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx);
            void generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::IndexType shotNumber);
            
            void getAcquisitionMat(Configuration::Configuration const &config, scai::lama::DenseMatrix<ValueType> &acqMat) const;

            scai::IndexType getNumShots();

          private:
            Seismogram<ValueType> signals; //!< Source signals
            suHandler<ValueType> su;

            /* Requiered acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_type; //!< Type of wavelet: 1==Synthetic

            /* Optional acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_shape;  //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            scai::lama::DenseVector<ValueType> wavelet_fc;     //!< Center frequency of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_amp;    //!< Amplitude of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_tshift; //!< Time shift of synthetic wavelet

            scai::IndexType numShots; //!< number of shots =1 for simultaneous execution =Number of Sources for serial execution
            bool wavelet_type_flag_2 = false; // flag if wavelet type 2 is used
            bool wavelet_type_flag_3 = false; // flag if wavelet type 3 is used

            void initOptionalAcquisitionParameter(scai::IndexType numParameter, scai::IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx) override;
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;

            void copySignalsToSeismogramHandler();

            void allocateSeismogram(scai::IndexType NT, scai::dmemo::DistributionPtr dist_traces, scai::hmemo::ContextPtr ctx);
            void generateSyntheticSignal(scai::IndexType SourceLocal, scai::IndexType NT, ValueType DT);
            void readSignalFromFile(Configuration::Configuration const &config, scai::IndexType SourceLocal, scai::IndexType numSourceRead);
        };
    }
}
