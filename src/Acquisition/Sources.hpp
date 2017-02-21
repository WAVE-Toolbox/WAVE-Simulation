#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>

#include "AcquisitionGeometry.hpp"
#include "SourceSignal/all.hpp"

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

            //! \brief Default destructor
            ~Sources(){};

            void init(scai::lama::DenseMatrix<ValueType> acquisition_temp,Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);

            void writeSignalsToFileRaw(std::string const &filename) const;

            void generateSignals(IndexType NT, ValueType DT, scai::hmemo::ContextPtr ctx);

          private:
            Seismogram<ValueType> signals; //!< Source signals

            /* Requiered acquisition Settings */
            scai::lama::DenseVector<IndexType> wavelet_type; //!< Type of wavelet: 1==Synthetic

            /* Optional acquisition Settings */
            scai::lama::DenseVector<IndexType> wavelet_shape;  //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            scai::lama::DenseVector<ValueType> wavelet_fc;     //!< Center frequency of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_amp;    //!< Amplitude of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_tshift; //!< Time shift of synthetic wavelet

            void initOptionalAcquisitionParameter(IndexType numParameter, IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx) override;
            void checkRequiredNumParameter(IndexType numParameterCheck) override;

            void copySignalsToSeismogramHandler();

            void allocateSeismogram(IndexType NT, scai::dmemo::DistributionPtr dist_traces, scai::hmemo::ContextPtr ctx);
            void generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT);
        };
    }
}
