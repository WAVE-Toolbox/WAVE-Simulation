#include "Sources.hpp"

using namespace scai;

/*! \brief Constructor based on the configuration class and the distribution of the wavefields. This constructor will call an initialisation function
 *
 * If runSimultaneousShots is true all sources are initialized. 
 * 
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
KITGPI::Acquisition::Sources<ValueType>::Sources(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    if (config.get<bool>("runSimultaneousShots")) {
        init(config, ctx, dist_wavefield);
        numShots = 1;
    } else {
        scai::lama::DenseMatrix<ValueType> acquisition_temp;
        acquisition_temp.readFromFile(config.get<std::string>("SourceFilename"));
        numShots = acquisition_temp.getNumRows();
    }
}

/*! \brief Init based on the configuration class and the distribution of the wavefields. This Init will read the acquistion of a single source from the Sourcefile
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param shotNumber Number of the source in the Source File
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, IndexType shotNumber)
{
    /* Read shotNumber row ofacquisition matrix */
    scai::lama::DenseStorage<ValueType> test;
    test.readFromFile(config.get<std::string>("SourceFilename"), shotNumber, 1);
    scai::lama::DenseMatrix<ValueType> acquisition_temp(test);

    this->init(acquisition_temp, config, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields. This function will read the acquistion from the Sourcefile
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    /* Read acquisition matrix */
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(config.get<std::string>("SourceFilename"));

    this->init(acquisition_temp, config, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 * acquistion matrix:
 * |           | X | Y | Z | SOURCE_TYPE | WAVELET_TYPE | WAVELET_SHAPE | FC | AMP | TShift |
 * | ----------| --| --| --| ----------- | -------------| --------------| ---| ----| ------ |
 * | 1. source |   |   |   |             |              |               |    |     |        |
 *
 * 
 \param config Configuration class, which is used to derive all requiered parameters
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 \param acquisition_matrix Dense Matrix which holds number of sources rows and number of source parameters columns
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    IndexType NT = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    /* Read acquisition from file */
    this->setAcquisition(acquisition_matrix, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), dist_wavefield, ctx);

    /* init seismogram handler */
    this->initSeismogramHandler(NT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setNormalizeTraces(config.get<IndexType>("NormalizeTraces"));

    /* Generate Signals */
    generateSignals(NT, config.get<ValueType>("DT"), ctx);
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

/*! \brief Generation of the source signals
 *
 * Allocation and calculation of the source signals accordingly to the source parameter vectors.
 * The calculation is performed locally on each node.
 *
 \param NT Number of time steps
 \param DT Time step interval
 \param ctx context
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::generateSignals(IndexType NT, ValueType DT, scai::hmemo::ContextPtr ctx)
{
    SCAI_ASSERT(this->getNumParameter() >= 5, "Number of source parameters < 5. Cannot generate signals. ");
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

        default:
            COMMON_THROWEXCEPTION("Unkown wavelet type ")
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

    SCAI_ASSERT(this->getNumParameter() >= 9, "Number of source parameters <= 9. Cannot generate synthetic signals. ");

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
    if (numParameter > 5) {
        wavelet_shape.allocate(numTracesGlobal);
        wavelet_fc.allocate(numTracesGlobal);
        wavelet_amp.allocate(numTracesGlobal);
        wavelet_tshift.allocate(numTracesGlobal);
    }

    lama::DenseVector<ValueType> tmp;   // needed for conversion

    /* Save source configurations from acquisition matrix in vectors */
    acquisition.getRow(tmp, 4);
    wavelet_type = lama::cast<IndexType>( tmp );
    if (numParameter > 5) {
        acquisition.getRow(tmp, 5);
        wavelet_shape = lama::cast<IndexType>( tmp );
        acquisition.getRow(wavelet_fc, 6);
        acquisition.getRow(wavelet_amp, 7);
        acquisition.getRow(wavelet_tshift, 8);
    }

    /* Redistribute source parameter vectors to corresponding processes */
    wavelet_type.redistribute(dist_wavefield_traces);
    wavelet_type.setContextPtr(ctx);
    if (numParameter > 5) {
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

/*! \brief Write source signals to file
 *
 \param filename Filename to write source signals
 */
template <typename ValueType>
void KITGPI::Acquisition::Sources<ValueType>::writeSignalsToFileRaw(std::string const &filename) const
{
    signals.writeToFileRaw(filename);
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

    /* Signals matix is row distributed according to dist_wavefield_sources, No column distribution */
    signals.allocate(ctx, dist_traces, NT);
    signals.setCoordinates(this->getCoordinates());
    signals.setContextPtr(ctx);
}

template class KITGPI::Acquisition::Sources<double>;
template class KITGPI::Acquisition::Sources<float>;
