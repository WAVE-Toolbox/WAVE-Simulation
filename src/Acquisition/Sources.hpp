#pragma once

#include <scai/dmemo.hpp>
#include <scai/lama.hpp>
#include <vector>
#include <scai/common/Walltime.hpp>

#include "AcquisitionGeometry.hpp"
#include "AcquisitionSettings.hpp"
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
            //! \brief Default destructor
            ~Sources(){};
            void init(std::vector<sourceSettings<ValueType>> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield);
            void init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseMatrix<ValueType> &signalMatrix);

            void generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, std::vector<scai::IndexType> rowinds);

            void getAcquisitionSettings(Configuration::Configuration const &config, ValueType shotIncr);
            void calcSourceSettingsEncode(scai::dmemo::CommunicatorPtr commAll, Configuration::Configuration const &config);
            void calcUniqueShotInds(scai::dmemo::CommunicatorPtr commAll, Configuration::Configuration const &config, std::vector<IndexType> &shotHistory, IndexType maxcount);
            void writeShotIndIncr(Configuration::Configuration const &config, std::vector<IndexType> uniqueShotNos);
            void writeSourceEncode(Configuration::Configuration const &config, std::string filename);
            scai::lama::DenseMatrix<ValueType> getsourcesignal();
            void setsourcesignal(scai::lama::DenseMatrix<ValueType> setsourcesignal);
            std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> getSourceSettings();
            std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> getSourceSettingsEncode();
            std::vector<IndexType> getShotIndIncr();
            std::vector<IndexType> getUniqueShotInds();

            using AcquisitionGeometry<ValueType>::isSeismic;
            
          private:
            Seismogram<ValueType> signals; //!< Source signals
            suHandler<ValueType> su;

            /* Requiered acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_type; //!< Type of wavelet: 1==Synthetic

            /* Optional acquisition Settings */
            scai::lama::DenseVector<scai::IndexType> wavelet_shape; //!< Shape of wavelet: 1==Ricker,2==Sinw,3==sin^3,4==FGaussian,5==Spike,6==integral sin^3
            scai::lama::DenseVector<ValueType> wavelet_fc;          //!< Center frequency of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_amp;         //!< Amplitude of synthetic wavelet
            scai::lama::DenseVector<ValueType> wavelet_tshift;      //!< Time shift of synthetic wavelet

            std::vector<sourceSettings<ValueType>> sourceSettingsShotIncr; // sourceSettings of selected shots
            std::vector<sourceSettings<ValueType>> sourceSettingsEncode; // sourceSettings of the encoded shots
            std::vector<IndexType> shotIndIncr;   // shot indices selected by shot increment
            std::vector<IndexType> uniqueShotInds;  // shot indices of random shots
    
            bool wavelet_type_flag_2 = false; // flag if wavelet type 2 is used
            bool wavelet_type_flag_3 = false; // flag if wavelet type 3 is used

            void initOptionalAcquisitionParameter(scai::IndexType numParameter, scai::IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx) override;
            void initOptionalAcquisitionParameter(std::vector<sourceSettings<ValueType>> allSettings, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx);
            void checkRequiredNumParameter(scai::IndexType numParameterCheck) override;
            void acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<sourceSettings<ValueType>> &allSettings, scai::dmemo::DistributionPtr dist_wavefield);

            void copySignalsToSeismogramHandler(std::vector<scai::IndexType> sourceNos);

            void allocateSeismogram(scai::IndexType NT, scai::dmemo::DistributionPtr dist_traces, scai::hmemo::ContextPtr ctx);
            void generateSyntheticSignal(scai::IndexType SourceLocal, scai::IndexType NT, ValueType DT);
            void readSignalFromFile(Configuration::Configuration const &config, scai::IndexType SourceLocal, scai::IndexType rowNumber);
        };
    }
}
