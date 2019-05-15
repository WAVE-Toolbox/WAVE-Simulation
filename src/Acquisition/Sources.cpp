#include "Sources.hpp"

using namespace scai;

/*! \brief Constructor based on the configuration class and the distribution of the wavefields. This constructor will call an initialisation function
 *
 * If runSimultaneousShots is true all sources are initialized. 
 * 
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
KITGPI::Acquisition::Sources<ValueType>::Sources(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    init(config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields. This Init will read the acquistion of a single source from the Sourcefile or from SU. If the acquisition should be read from SU all sources have to be initialized from SU first.
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param shotNumber Number of the source in the Source File
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, IndexType shotNumber)
{
    /* Read shotNumber row of acquisition matrix */
    std::vector<sourceSettings<ValueType>> allSettings;
    sourceSettings<ValueType> settings;
    
    if (!config.get<bool>("initSourcesFromSU")) {
        readSettings(settings, config.get<std::string>("SourceFilename") + ".txt", shotNumber);
    } else {
        settings = su.getSourceSettings(shotNumber);
    }

    allSettings.push_back(settings);
    this->init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init of all sources based on the configuration class and the distribution of the wavefields. This function will read the acquistion from the Sourcefile or from SU.
 * If SU is selected the following header words have to be set: sx, sy, sdepth, scalco and scalel.
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    /* Read acquisition matrix */
    std::vector<sourceSettings<ValueType>> allSettings;

    if (config.get<bool>("initSourcesFromSU")) {
        su.buildAcqMatrixSource(config.get<std::string>("SourceSignalFilename"), modelCoordinates.getDH());
        allSettings = su.getSourceSettingsVec();
    } else
        readAllSettings(allSettings, config.get<std::string>("SourceFilename") + ".txt");

    this->init(allSettings, config, modelCoordinates, ctx, dist_wavefield);

    if (config.get<bool>("runSimultaneousShots")) {
        numShots = 1;
    } else {
        numShots = allSettings.size();
    }
}


/*! \brief Init of all shots based on the configuration class and the distribution of the wavefields
 *
 * acquistion matrix:
 * |           | X | Y | Z | SOURCE_TYPE | WAVELET_TYPE | WAVELET_SHAPE | FC | AMP | TShift |
 * | ----------| --| --| --| ----------- | -------------| --------------| ---| ----| ------ |
 * | 1. source |   |   |   |             |              |               |    |     |        |
 *
 * 
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param acquisition_matrix Dense Matrix which holds number of sources rows and number of source parameters columns
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(std::vector<sourceSettings<ValueType>> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    /*reset seismograms. This is necessary when init will be called multiple times*/
    this->getSeismogramHandler().resetSeismograms();

    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* Read acquisition from file */
    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);
    initOptionalAcquisitionParameter(allSettings, dist_wavefield, ctx);
    
    /* init seismogram handler */
    this->initSeismogramHandler(NT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<IndexType>("NormalizeTraces"));

    /* Generate Signals */
    generateSignals(config, ctx);
    copySignalsToSeismogramHandler();
}

/*! \brief Init of a single shot based on the configuration class and the distribution of the wavefields
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param acquisition_matrix Dense Matrix which holds number of sources rows and number of source parameters columns
 \param shotNumber Shot number
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(std::vector<sourceSettings<ValueType>> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber)
{
    /*reset seismograms. This is necessary when init will be called multiple times*/
    this->getSeismogramHandler().resetSeismograms();

    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* Read acquisition from file */
    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);
    initOptionalAcquisitionParameter(allSettings, dist_wavefield, ctx);
    
    /* init seismogram handler */
    this->initSeismogramHandler(NT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<IndexType>("NormalizeTraces"));

    /* Generate Signals */
    generateSignals(config, ctx, shotNumber);
    copySignalsToSeismogramHandler();
}

/*! \brief Init with a signal matrix
 \param config Configuration class, which is used to derive all requiered parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param acquisition_matrix Dense Matrix which holds number of sources rows and number of source parameters columns
 \param signalMatrix Signal matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseMatrix<ValueType> &signalMatrix)
{
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
    this->getSeismogramHandler().setNormalizeTraces(config.get<IndexType>("NormalizeTraces"));

    /* Generate Signals */
    allocateSeismogram(NT, this->getSeismogramTypes().getDistributionPtr(), ctx);
    signals.setDT(DT);
    signalMatrix.redistribute(signals.getData().getRowDistributionPtr(), signals.getData().getColDistributionPtr());
    signals.getData() = signalMatrix;
    copySignalsToSeismogramHandler();
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Sources<ValueType>::getNumShots()
{
    return numShots;
}

/*! \brief Generation of the source signals for all shots
 *
 * Allocation and calculation of the source signals accordingly to the source parameter vectors.
 * The calculation is performed locally on each node.
 *
 \param NT Number of time steps
 \param DT Time step interval
 \param ctx context
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx)
{
    ValueType DT = config.get<ValueType>("DT");
    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / DT) + 0.5);

    SCAI_ASSERT_GT_DEBUG(NT, 0, "NT must be positive");
    SCAI_ASSERT_GT_DEBUG(DT, 0, "DT must be positive");

    allocateSeismogram(NT, this->getSeismogramTypes().getDistributionPtr(), ctx);

    signals.setDT(DT);

    auto read_wavelet_type_LA = hmemo::hostReadAccess(wavelet_type.getLocalValues());

    IndexType wavelet_type_i;

    for (IndexType i = 0; i < this->getNumTracesLocal(); i++) {

        /* Cast to IndexType */
        wavelet_type_i = read_wavelet_type_LA[i];

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
            readSignalFromFile(config, i, i);
            break;

        default:
            COMMON_THROWEXCEPTION("Unkown wavelet type ")
            break;
        }
    }
}

/*! \brief Generation of the source signal for a single shot
 *
 * Allocation and calculation of the source signals accordingly to the source parameter vectors.
 * The calculation is performed locally on each node.
 *
 \param NT Number of time steps
 \param DT Time step interval
 \param ctx context
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::IndexType shotNumber)
{
    ValueType DT = config.get<ValueType>("DT");
    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / DT) + 0.5);

    SCAI_ASSERT_GT_DEBUG(NT, 0, "NT must be positive");
    SCAI_ASSERT_GT_DEBUG(DT, 0, "DT must be positive");

    allocateSeismogram(NT, this->getSeismogramTypes().getDistributionPtr(), ctx);

    signals.setDT(DT);

    IndexType wavelet_type_i = wavelet_type[0];

    switch (wavelet_type_i) {
    case 1:
        /* Synthetic wavelet */
        generateSyntheticSignal(0, NT, DT);
        break;
    case 2: // read one source signal for all sources
        wavelet_type_flag_2 = true;
        SCAI_ASSERT(!wavelet_type_flag_3, "Combination of wavelet type 2 and 3 not supported");
        readSignalFromFile(config, 0, 0);
        break;
    case 3:
        wavelet_type_flag_3 = true; // read a source signal for each source
        SCAI_ASSERT(!wavelet_type_flag_2, "Combination of wavelet type 2 and 3 not supported");
        readSignalFromFile(config, 0, shotNumber);
        break;

    default:
        COMMON_THROWEXCEPTION("Unkown wavelet type ")
        break;
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
        SourceSignal::Spike<ValueType>(signalVector, NT, DT, 0, wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    case 6:
        /* integral sin3 signal */
        SourceSignal::IntgSinThree<ValueType>(signalVector, NT, DT, wavelet_fc.getLocalValues()[SourceLocal], wavelet_amp.getLocalValues()[SourceLocal], wavelet_tshift.getLocalValues()[SourceLocal]);
        break;

    default:
        COMMON_THROWEXCEPTION("Unkown wavelet shape ")
        break;
    }

    hmemo::HArray<ValueType> localsignal = signalVector.getLocalValues();

    lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
    signalsMatrix.setLocalRow(localsignal, SourceLocal, scai::common::BinaryOp::COPY);
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::readSignalFromFile(Configuration::Configuration const &config, scai::IndexType SourceLocal, scai::IndexType numSourceRead)
{
    std::string signalFilename = config.get<std::string>("SourceSignalFilename");

    if (config.get<bool>("initSourcesFromSU")) {
        scai::lama::DenseVector<ValueType> singleSignal;
        scai::IndexType numSourceReadLocal;

        su.locateTrace(signalFilename, numSourceReadLocal, numSourceRead);
        suHandler<ValueType>::readSingleDataSU(signalFilename, singleSignal, numSourceReadLocal);

        SCAI_ASSERT(singleSignal.size() == signals.getData().getNumColumns(), "Source signal has invalid length");

        lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
        signalsMatrix.setRow(singleSignal, SourceLocal, scai::common::BinaryOp::COPY);
    } else {
        signalFilename += ".mtx";
        scai::lama::DenseStorage<ValueType> signalsMatrixTmp;
        signalsMatrixTmp.readFromFile(signalFilename, numSourceRead, 1);

        SCAI_ASSERT(signalsMatrixTmp.getNumColumns() == signals.getData().getNumColumns(), "Source signal has invalid length");

        hmemo::HArray<ValueType> localsignal = signalsMatrixTmp.getValues();
        lama::DenseMatrix<ValueType> &signalsMatrix = signals.getData();
        signalsMatrix.setLocalRow(localsignal, SourceLocal, scai::common::BinaryOp::COPY);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::checkRequiredNumParameter(IndexType numParameterCheck)
{

    if (numParameterCheck < 5 || numParameterCheck > 9) {
        COMMON_THROWEXCEPTION("Source acquisition file has an unkown format ")
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

        seismograms.getSeismogram(static_cast<SeismogramType>(tempIndexType)).getData().setRow(temp, count[tempIndexType], scai::common::BinaryOp::COPY);

        ++count[tempIndexType];
    }

    SCAI_ASSERT_DEBUG(count[0] == seismograms.getNumTracesGlobal(SeismogramType::P), " Size mismatch ");
    SCAI_ASSERT_DEBUG(count[1] == seismograms.getNumTracesGlobal(SeismogramType::VX), " Size mismatch ");
    SCAI_ASSERT_DEBUG(count[2] == seismograms.getNumTracesGlobal(SeismogramType::VY), " Size mismatch ");
    SCAI_ASSERT_DEBUG(count[3] == seismograms.getNumTracesGlobal(SeismogramType::VZ), " Size mismatch ");
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
 * Uses configuration to determine if sources are initialized by MTX or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::getAcquisitionMat(Configuration::Configuration const &config, std::vector<sourceSettings<ValueType>> &allSourceSettings)
{
    if (config.get<bool>("initSourcesFromSU"))
        allSourceSettings = su.getSourceSettingsVec();
    else
        readAllSettings(allSourceSettings, config.get<std::string>("SourceFilename") + ".mtx");
}

template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<sourceSettings<ValueType>> &allSettings, scai::dmemo::DistributionPtr dist_wavefield) {
    scai::IndexType nRows = acqMat.getNumRows();
    allSettings.reserve(nRows);
    if (dist_wavefield->getCommunicator().getRank() == 0) {
        auto read_acquisition_temp_HA = hostReadAccess(acqMat.getLocalStorage().getData());
        for (IndexType row = 0; row < nRows; row++) {
            allSettings[row].sourceCoords.x = read_acquisition_temp_HA[row*9+0];
            allSettings[row].sourceCoords.y = read_acquisition_temp_HA[row*9+1];
            allSettings[row].sourceCoords.z = read_acquisition_temp_HA[row*9+2];
            allSettings[row].sourceType = read_acquisition_temp_HA[row*9+3];
            allSettings[row].waveletType = read_acquisition_temp_HA[row*9+4];
            allSettings[row].waveletShape = read_acquisition_temp_HA[row*9+5];
            allSettings[row].fc = read_acquisition_temp_HA[row*9+6];
            allSettings[row].amp = read_acquisition_temp_HA[row*9+7];
            allSettings[row].tShift = read_acquisition_temp_HA[row*9+8];
        }
    }
}


template class KITGPI::Acquisition::Sources<double>;
template class KITGPI::Acquisition::Sources<float>;
