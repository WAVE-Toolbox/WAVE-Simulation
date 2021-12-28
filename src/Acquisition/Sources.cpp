
#include "Sources.hpp"
#include "../IO/IO.hpp"
#include "../IO/SUIO.hpp"

using namespace scai;

/*! \brief Init of a single shot based on the configuration class and the distribution of the wavefields
 \param allSettings vector of sourceSettings structs with settings for all shots
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(std::vector<sourceSettings<ValueType>> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{    
    /* identify seismic and EM wave */
    std::string equationType = config.get<std::string>("equationType");
    isSeismic = Common::checkEquationType<ValueType>(equationType);
    this->getSeismogramHandler().setIsSeismic(isSeismic);
    
    /*reset seismograms. This is necessary when init will be called multiple times*/
    this->getSeismogramHandler().resetSeismograms();

    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* Read acquisition from file */
    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);
    initOptionalAcquisitionParameter(allSettings, dist_wavefield, ctx);

    /* init seismogram handler */
    this->initSeismogramHandler(NT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));

    /* Generate Signals */
    std::vector<scai::IndexType> readrows;
    for (unsigned long i = 0; i < allSettings.size(); i++) {
        readrows.push_back(allSettings[i].row);
    }

    generateSignals(config, ctx, readrows);
    copySignalsToSeismogramHandler();
}

/*! \brief Init with a signal matrix
 \param acquisition_matrix Dense Matrix which holds number of sources rows and number of source parameters columns
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param signalMatrix Signal matrix
 */

// this function can be removed acq matrix is no longer supported maybe someone still uses it?

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseMatrix<ValueType> &signalMatrix)
{
    /* identify seismic and EM wave */
    std::string equationType = config.get<std::string>("equationType");
    isSeismic = Common::checkEquationType<ValueType>(equationType);
    this->getSeismogramHandler().setIsSeismic(isSeismic);
    
    ValueType DT = config.get<ValueType>("DT");
    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / DT) + 0.5);

    /* Read acquisition from file */
    std::vector<sourceSettings<ValueType>> allSettings;
    acqMat2settings(acquisition_matrix, allSettings, dist_wavefield);
    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);
    initOptionalAcquisitionParameter(allSettings, dist_wavefield, ctx);

    /* init seismogram handler */
    this->initSeismogramHandler(NT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));

    /* Generate Signals */
    allocateSeismogram(NT, this->getSeismogramTypes().getDistributionPtr(), ctx);
    signals.setDT(DT);
    signalMatrix.redistribute(signals.getData().getRowDistributionPtr(), signals.getData().getColDistributionPtr());
    signals.getData() = signalMatrix;
    copySignalsToSeismogramHandler();
}

/*! \brief Generation of the source signal for a single shot
 *
 * Allocation and calculation of the source signals accordingly to the source parameter vectors.
 * The calculation is performed locally on each node.
 *
 \param config Configuration file
 \param ctx context
 \param rowinds vector with rows to read from source signal file corresponding to lines in sources.txt
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, std::vector<scai::IndexType> rowinds)
{
    ValueType DT = config.get<ValueType>("DT");
    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / DT) + 0.5);

    SCAI_ASSERT_GT_DEBUG(NT, 0, "NT must be positive");
    SCAI_ASSERT_GT_DEBUG(DT, 0, "DT must be positive");

    allocateSeismogram(NT, this->getSeismogramTypes().getDistributionPtr(), ctx);

    signals.setDT(DT);

    auto read_wavelet_type_LA = hmemo::hostReadAccess(wavelet_type.getLocalValues());

    IndexType wavelet_type_i;

    scai::dmemo::DistributionPtr shotdist = this->getSeismogramTypes().getDistributionPtr();

    for (IndexType i = 0; i < this->getNumTracesLocal(); i++) {
        /* Cast to IndexType */
        wavelet_type_i = read_wavelet_type_LA[i];

        // get global indices for local ones
        IndexType glob_i = shotdist->local2Global(i);

        switch (wavelet_type_i) {
        case 1:
            /* Synthetic wavelet */
            generateSyntheticSignal(i, NT, DT);
            break;
        case 2: // read one source signal for all sources
            wavelet_type_flag_2 = true;
            SCAI_ASSERT(!wavelet_type_flag_3, "Combination of wavelet type 2 and 3 not supported");
            readSignalFromFile(config, i, 0);
            break;
        case 3:
            wavelet_type_flag_3 = true; // read a source signal for each source
            SCAI_ASSERT(!wavelet_type_flag_2, "Combination of wavelet type 2 and 3 not supported");
            readSignalFromFile(config, i, rowinds[glob_i]);
            break;

        default:
            COMMON_THROWEXCEPTION("Unknown wavelet type ")
            break;
        }
    }
}

/*! \brief Generation of synthetic source signals
 *
 * Calculation of a synthetic source signal accordingly to the source parameter vectors for the given local source number.
 * Uses the entries of the wavelet_shape vector to determine the shape of the wavelet.
 *
 \param SourceLocal Number of the local source
 \param NT Number of time steps
 \param DT Time step interval
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSyntheticSignal(IndexType SourceLocal, IndexType NT, ValueType DT)
{
    lama::DenseVector<ValueType> signalVector;
    signalVector.allocate(NT);

    /* Cast to IndexType */
    IndexType wavelet_shape_i = wavelet_shape.getLocalValues()[SourceLocal];

    switch (wavelet_shape_i) {
    case 1:
        /* Ricker */
        SourceSignal::Ricker<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 2:
        /* combination of sin signals */
        SourceSignal::SinW<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 3:
        /* sin3 signal */
        SourceSignal::SinThree<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 4:
        /* First derivative of a Gaussian (FGaussian) */
        SourceSignal::FGaussian<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 5:
        /* Spike signal */
        SourceSignal::Spike<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 6:
        /* integral sin3 signal */
        SourceSignal::IntgSinThree<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;
        
    case 7:
        /* Ricker_GprMax */
        SourceSignal::Ricker_GprMax<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;
        
    case 8:
        /* Berlage */
        SourceSignal::Berlage<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;
        
    default:
        COMMON_THROWEXCEPTION("Unknown wavelet shape ")
        break;
    }

    hmemo::HArray<ValueType> localsignal = signalVector.getLocalValues();

    lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
    signalsMatrix.setLocalRow(localsignal, SourceLocal, scai::common::BinaryOp::COPY);
}

/*! \brief Read source signal from file
 *
 * Read source signal from file
 *
 \param config Configuration file
 \param SourceLocal Number of the local source
 \param rowNumber Number of the row to read in
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::readSignalFromFile(Configuration::Configuration const &config, scai::IndexType SourceLocal, scai::IndexType rowNumber)
{
    std::string signalFilename = config.get<std::string>("SourceSignalFilename");

    if (config.get<bool>("initSourcesFromSU")) {
        scai::lama::DenseVector<ValueType> singleSignal;

        SUIO::readSingleDataSU(signalFilename, singleSignal, rowNumber);

        SCAI_ASSERT(singleSignal.size() == signals.getData().getNumColumns(), "Source signal has invalid length");

        hmemo::HArray<ValueType> localsignal = singleSignal.getLocalValues();
        lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
        signalsMatrix.setLocalRow(localsignal, SourceLocal, scai::common::BinaryOp::COPY);
    } else {
        hmemo::HArray<ValueType> localsignal = IO::readMatrix<ValueType>(signalFilename, rowNumber, config.get<IndexType>("SeismogramFormat"));
        lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
        signalsMatrix.setLocalRow(localsignal, SourceLocal, scai::common::BinaryOp::COPY);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::checkRequiredNumParameter(IndexType numParameterCheck)
{
    if (numParameterCheck < 6 || numParameterCheck > 10) {
        COMMON_THROWEXCEPTION("Source acquisition file has an unknown format ")
    }
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::initOptionalAcquisitionParameter(IndexType numParameter, IndexType numTracesGlobal, scai::lama::DenseMatrix<ValueType> acquisition, scai::dmemo::DistributionPtr dist_wavefield_traces, scai::hmemo::ContextPtr ctx)
{
    /* Allocate source parameter vectors on all processes */
    wavelet_type.allocate(numTracesGlobal);
    if (numParameter > 4) {
        wavelet_shape.allocate(numTracesGlobal);
        wavelet_fc.allocate(numTracesGlobal);
        wavelet_amp.allocate(numTracesGlobal);
        wavelet_tshift.allocate(numTracesGlobal);
    }

    lama::DenseVector<ValueType> tmp; // needed for conversion

    /* Save source configurations from acquisition matrix in vectors */
    acquisition.getRow(tmp, 4);
    wavelet_type = lama::cast<IndexType>(tmp);
    if (numParameter > 4) {
        acquisition.getRow(tmp, 5);
        wavelet_shape = lama::cast<IndexType>(tmp);
        acquisition.getRow(wavelet_fc, 6);
        acquisition.getRow(wavelet_amp, 7);
        acquisition.getRow(wavelet_tshift, 8);
    }

    /* Redistribute source parameter vectors to corresponding processes */
    wavelet_type.redistribute(dist_wavefield_traces);
    wavelet_type.setContextPtr(ctx);
    if (numParameter > 4) {
        wavelet_shape.redistribute(dist_wavefield_traces);
        wavelet_fc.redistribute(dist_wavefield_traces);
        wavelet_amp.redistribute(dist_wavefield_traces);
        wavelet_tshift.redistribute(dist_wavefield_traces);

        wavelet_shape.setContextPtr(ctx);
        wavelet_fc.setContextPtr(ctx);
        wavelet_amp.setContextPtr(ctx);
        wavelet_tshift.setContextPtr(ctx);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::initOptionalAcquisitionParameter(std::vector<sourceSettings<ValueType>> allSettings, scai::dmemo::DistributionPtr dist_wavefield, scai::hmemo::ContextPtr ctx)
{
    /* Allocate source parameter vectors on all processes */
    scai::IndexType numTracesGlobal = this->getNumTracesGlobal();
    scai::dmemo::DistributionPtr dist_master_numTracesGlobal(new scai::dmemo::CyclicDistribution(numTracesGlobal, numTracesGlobal, dist_wavefield->getCommunicatorPtr()));
    scai::dmemo::DistributionPtr dist_wavefield_traces = this->getSeismogramTypes().getDistributionPtr();

    wavelet_type.allocate(dist_master_numTracesGlobal);
    wavelet_shape.allocate(dist_master_numTracesGlobal);
    wavelet_fc.allocate(dist_master_numTracesGlobal);
    wavelet_amp.allocate(dist_master_numTracesGlobal);
    wavelet_tshift.allocate(dist_master_numTracesGlobal);

    if (dist_wavefield->getCommunicator().getRank() == 0) {
        /* Get writeAccess to coordinates vector (local) */
        auto write_waveletType_LA = hostWriteAccess(wavelet_type.getLocalValues());
        auto write_waveletShape_LA = hostWriteAccess(wavelet_shape.getLocalValues());
        auto write_waveletFc_LA = hostWriteAccess(wavelet_fc.getLocalValues());
        auto write_amp_LA = hostWriteAccess(wavelet_amp.getLocalValues());
        auto write_tshift_LA = hostWriteAccess(wavelet_tshift.getLocalValues());

        /* 2. Calculate 1-D coordinates from 3-D coordinates */
        for (scai::IndexType i = 0; i < numTracesGlobal; i++) {
            write_waveletType_LA[i] = allSettings[i].waveletType;
            write_waveletShape_LA[i] = allSettings[i].waveletShape;
            write_waveletFc_LA[i] = allSettings[i].fc;
            write_amp_LA[i] = allSettings[i].amp;
            write_tshift_LA[i] = allSettings[i].tShift;
        }

        write_waveletType_LA.release();
        write_waveletShape_LA.release();
        write_waveletFc_LA.release();
        write_amp_LA.release();
        write_tshift_LA.release();
    }

    /* Redistribute source parameter vectors to corresponding processes */
    wavelet_type.redistribute(dist_wavefield_traces);
    wavelet_shape.redistribute(dist_wavefield_traces);
    wavelet_fc.redistribute(dist_wavefield_traces);
    wavelet_amp.redistribute(dist_wavefield_traces);
    wavelet_tshift.redistribute(dist_wavefield_traces);

    wavelet_type.setContextPtr(ctx);
    wavelet_shape.setContextPtr(ctx);
    wavelet_fc.setContextPtr(ctx);
    wavelet_amp.setContextPtr(ctx);
    wavelet_tshift.setContextPtr(ctx);
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::copySignalsToSeismogramHandler()
{
    IndexType tempIndexType;
    lama::DenseVector<ValueType> temp;
    SeismogramHandler<ValueType> &seismograms = this->getSeismogramHandler();
    IndexType count[NUM_ELEMENTS_SEISMOGRAMTYPE] = {0, 0, 0, 0};

    /* Copy data to the seismogram handler */
    for (IndexType i = 0; i < this->getNumTracesGlobal(); ++i) {
        tempIndexType = this->getSeismogramTypes().getValue(i) - 1;

        signals.getData().getRow(temp, i);

        if (isSeismic)
            seismograms.getSeismogram(static_cast<SeismogramType>(tempIndexType)).getData().setRow(temp, count[tempIndexType], scai::common::BinaryOp::COPY);
        else
            seismograms.getSeismogram(static_cast<SeismogramTypeEM>(tempIndexType)).getData().setRow(temp, count[tempIndexType], scai::common::BinaryOp::COPY);

        ++count[tempIndexType];
    }
    if (isSeismic) {
        SCAI_ASSERT_DEBUG(count[0] == seismograms.getNumTracesGlobal(SeismogramType::P), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[1] == seismograms.getNumTracesGlobal(SeismogramType::VX), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[2] == seismograms.getNumTracesGlobal(SeismogramType::VY), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[3] == seismograms.getNumTracesGlobal(SeismogramType::VZ), " Size mismatch ");
    } else {
        SCAI_ASSERT_DEBUG(count[0] == seismograms.getNumTracesGlobal(SeismogramTypeEM::EZ), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[1] == seismograms.getNumTracesGlobal(SeismogramTypeEM::EX), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[2] == seismograms.getNumTracesGlobal(SeismogramTypeEM::EY), " Size mismatch ");
        SCAI_ASSERT_DEBUG(count[3] == seismograms.getNumTracesGlobal(SeismogramTypeEM::HZ), " Size mismatch ");
    }
}

/*! \brief get source signal
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType> KITGPI::Acquisition::Sources<ValueType>::getsourcesignal()
{
    lama::DenseMatrix<ValueType> signal_out = signals.getData();
    return (signal_out);
}

/*! \brief set source signal
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::setsourcesignal(lama::DenseMatrix<ValueType> setsourcesignal)
{
    signals.getData() = setsourcesignal;
    copySignalsToSeismogramHandler();
}

/*! \brief Allocation of the source signals matrix
 *
 * Allocation of the source signals matrix based on an already defined source distribution and the number of time steps.
 * The source signal matrix is allocated based on the distributions.
 *
 \param NT Number of time steps
 \param ctx context
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::allocateSeismogram(IndexType NT, scai::dmemo::DistributionPtr dist_traces, scai::hmemo::ContextPtr ctx)
{
    SCAI_ASSERT_DEBUG(NT > 0, "NT<=0");
    if (dist_traces == NULL) {
        COMMON_THROWEXCEPTION("Row distribution of sources (dist_wavefield_sources) is not set!")
    }

    /* Signals matrix is row distributed according to dist_wavefield_sources, No column distribution */
    signals.allocate(ctx, dist_traces, NT);
    signals.setCoordinates(this->get1DCoordinates());
    signals.setContextPtr(ctx);
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param allSettings sourceSettings
 \param shotIncr shot increments in meters
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<sourceSettings<ValueType>> &allSettings, std::vector<IndexType> &shotIndIncr, ValueType shotIncr)
{
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings; 
    if (config.get<bool>("initSourcesFromSU"))
        su.readAllSettingsFromSU(sourceSettings, config.get<std::string>("SourceFilename"), config.get<ValueType>("DH"));
    else
        readAllSettings(sourceSettings, config.get<std::string>("SourceFilename") + ".txt");
    IndexType numshots = sourceSettings.size();
    shotIndIncr.clear();
    if (numshots > 1 && shotIncr > config.get<ValueType>("DH")) {
        std::vector<IndexType> sourceLength; 
        if (abs(sourceSettings[0].sourceCoords.x-sourceSettings[1].sourceCoords.x) > abs(sourceSettings[0].sourceCoords.y-sourceSettings[1].sourceCoords.y)) {
            for (IndexType shotInd = 0; shotInd < numshots; shotInd++) {
                sourceLength.push_back(sourceSettings[shotInd].sourceCoords.x);
            }
        } else if (abs(sourceSettings[0].sourceCoords.x-sourceSettings[1].sourceCoords.x) < abs(sourceSettings[0].sourceCoords.y-sourceSettings[1].sourceCoords.y)) {
            for (IndexType shotInd = 0; shotInd < numshots; shotInd++) {
                sourceLength.push_back(sourceSettings[shotInd].sourceCoords.y);
            }
        }
        IndexType numshotsIncr = 0;
        IndexType shotIncrInd = sourceLength[0];
        allSettings.clear();
        allSettings.push_back(sourceSettings[0]);
        shotIndIncr.push_back(0);
        IndexType shotIncrStep = floor(shotIncr / config.get<ValueType>("DH"));
        for (IndexType shotInd = 0; shotInd < numshots-1; shotInd++) {
            if (shotIncrInd >= sourceLength[shotInd] && shotIncrInd <= sourceLength[shotInd+1]) {
                if (shotIndIncr[numshotsIncr] != shotInd && abs(shotIncrInd - sourceLength[shotInd]) < abs(shotIncrInd - sourceLength[shotInd+1])) {
                    allSettings.push_back(sourceSettings[shotInd]);
                    shotIndIncr.push_back(shotInd);
                    numshotsIncr++; // to avoid repeat
                } else if (shotIndIncr[numshotsIncr] != shotInd+1 && abs(shotIncrInd - sourceLength[shotInd]) >= abs(shotIncrInd - sourceLength[shotInd+1])) {
                    allSettings.push_back(sourceSettings[shotInd+1]);
                    shotIndIncr.push_back(shotInd+1);
                    numshotsIncr++; // to avoid repeat
                }
                shotIncrInd += shotIncrStep;
                shotInd--; // to avoid that some locations are skipped
            } else if (shotIndIncr[numshotsIncr] != shotInd+1 && shotInd == numshots - 2 && abs(sourceLength[shotInd+1] - sourceLength[shotInd]) > abs(shotIncrInd - sourceLength[shotInd+1])) {
                allSettings.push_back(sourceSettings[shotInd+1]);
                shotIndIncr.push_back(shotInd+1);
            }
        } 
    } else {
        allSettings = sourceSettings;
        for (IndexType shotInd = 0; shotInd < numshots; shotInd++) {
            shotIndIncr.push_back(shotInd);
        }
    }
}

/*! \brief Write the shot indices and shot numbers
 *
 \param config Configuration
 \param shotIndIncr shot indices
 \param uniqueShotNos shot numbers
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::writeShotIndIncr(Configuration::Configuration const &config, std::vector<IndexType> shotIndIncr, std::vector<IndexType> uniqueShotNos)
{
    ValueType shotIncr = config.getAndCatch("shotIncr", 0.0);
    if (shotIncr > config.get<ValueType>("DH")) {
        SCAI_ASSERT_ERROR(shotIndIncr.size() == uniqueShotNos.size(), "shotIndIncr.size() != uniqueShotNos.size()"); // check whether shotIncr has been applied successfully.
        std::string filename;
        bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
        if (useStreamConfig) {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            filename = configBig.get<std::string>("SourceFilename") + ".shotIncr.txt";
        } else {
            filename = config.get<std::string>("SourceFilename") + ".shotIncr.txt";
        }
        IndexType numshotsIncr = shotIndIncr.size();
        std::ofstream outputFile; 
        outputFile.open(filename);
        outputFile << "# Shot indices (shotIncr = " << shotIncr << " m, numshots = " << numshotsIncr << ")\n"; 
        outputFile << "# Shot index | shot number\n"; 
        for (int shotInd = 0; shotInd < numshotsIncr; shotInd++) { 
            outputFile << std::setw(12) << shotIndIncr[shotInd]+1 << std::setw(12) << uniqueShotNos[shotInd] << "\n";
        }
    }
}

/*! \brief Gets the Acquisition settings from Matrix
 *
 \param acqMat Matrix defined by user
 \param allSettings sourceSettings
 \param dist_wavefield dist_wavefield
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<sourceSettings<ValueType>> &allSettings, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::IndexType nRows = acqMat.getNumRows();
    allSettings.reserve(nRows);
    if (dist_wavefield->getCommunicator().getRank() == 0) {
        auto read_acquisition_temp_HA = hostReadAccess(acqMat.getLocalStorage().getData());
        for (IndexType row = 0; row < nRows; row++) {
            allSettings[row].sourceNo = read_acquisition_temp_HA[row * 10 + 0];
            allSettings[row].sourceCoords.x = read_acquisition_temp_HA[row * 10 + 1];
            allSettings[row].sourceCoords.y = read_acquisition_temp_HA[row * 10 + 2];
            allSettings[row].sourceCoords.z = read_acquisition_temp_HA[row * 10 + 3];
            allSettings[row].sourceType = read_acquisition_temp_HA[row * 10 + 4];
            allSettings[row].waveletType = read_acquisition_temp_HA[row * 10 + 5];
            allSettings[row].waveletShape = read_acquisition_temp_HA[row * 10 + 6];
            allSettings[row].fc = read_acquisition_temp_HA[row * 10 + 7];
            allSettings[row].amp = read_acquisition_temp_HA[row * 10 + 8];
            allSettings[row].tShift = read_acquisition_temp_HA[row * 10 + 9];
        }
    }
}

template class KITGPI::Acquisition::Sources<double>;
template class KITGPI::Acquisition::Sources<float>;
