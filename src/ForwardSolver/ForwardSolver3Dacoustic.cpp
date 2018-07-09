#include "ForwardSolver3Dacoustic.hpp"
using namespace scai;

/*! \brief Initialitation of the ForwardSolver
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param wavefield Wavefields for the modelling
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, scai::hmemo::ContextPtr ctx, ValueType /*DT*/)
{
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefVX().getDistributionPtr();
    ;

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, derivatives, dist, ctx);
    }

    /* Initialisation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
}

/*! \brief resets PML (use after each modelling!)
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::resetCPML()
{
    if (useConvPML) {
        ConvPML.resetCPML();
    }
}

/*! \brief Initialitation of the boundary conditions
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{

    useFreeSurface = config.get<IndexType>("FreeSurface");

    /* Prepare Free Surface */
    if (useFreeSurface==1) {
        FreeSurface.init(dist, derivatives, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<ValueType>("DH"));
    }

    /* Prepare Damping Boundary */
    if (config.get<IndexType>("DampingBoundary") == 1) {
        useDampingBoundary = true;
        DampingBoundary.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"), useFreeSurface);
    }

    if (config.get<IndexType>("DampingBoundary") == 2) {
        useConvPML = true;
        ConvPML.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<IndexType>("DH"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("NPower"), config.get<ValueType>("KMaxCPML"), config.get<ValueType>("CenterFrequencyCPML"), config.get<ValueType>("VMaxCPML"), useFreeSurface);
    }
}

/*! \brief Running the 3-D acoustic foward solver
 *
 * Start the 3-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param tStart Counter start in for loop over time steps
 \param tEnd Counter end  in for loop over time steps
 \param DT Temporal Sampling intervall in seconds
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType t)
{

    SCAI_REGION("timestep");

    /* Get references to required modelparameter */
    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector<ValueType> const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector<ValueType> const &inverseDensityAverageZ = model.getInverseDensityAverageZ();

    /* Get references to required wavefields */
    lama::Vector<ValueType> &vX = wavefield.getRefVX();
    lama::Vector<ValueType> &vY = wavefield.getRefVY();
    lama::Vector<ValueType> &vZ = wavefield.getRefVZ();
    lama::Vector<ValueType> &p = wavefield.getRefP();

    /* Get references to required derivatives matrixes */
    lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    lama::Matrix<ValueType> const &Dzf = derivatives.getDzf();
    lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    lama::Matrix<ValueType> const &Dzb = derivatives.getDzb();
    lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
    lama::Matrix<ValueType> const &DyfFreeSurface = derivatives.getDyfFreeSurface();

    SourceReceiverImpl::FDTD3Dacoustic<ValueType> SourceReceiver(sources, receiver, wavefield);

        /* update velocity */
        update = Dxf * p;
        if (useConvPML) {
            ConvPML.apply_p_x(update);
        }
        update *= inverseDensityAverageX;
        vX += update;

        if (useFreeSurface == 1) {
            /* Apply image method */
            update = DyfFreeSurface * p;
        } else {
            update = Dyf * p;
        }

        if (useConvPML) {
            ConvPML.apply_p_y(update);
        }
        update *= inverseDensityAverageY;
        vY += update;

        update = Dzf * p;
        if (useConvPML) {
            ConvPML.apply_p_z(update);
        }
        update *= inverseDensityAverageZ;
        vZ += update;

        /* pressure update */
        update = Dxb * vX;
        if (useConvPML) {
            ConvPML.apply_vxx(update);
        }

        update_temp = Dyb * vY;
        if (useConvPML) {
            ConvPML.apply_vyy(update_temp);
        }
        update += update_temp;

        update_temp = Dzb * vZ;
        if (useConvPML) {
            ConvPML.apply_vzz(update_temp);
        }
        update += update_temp;

        update *= pWaveModulus;
        p += update;

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(p, vX, vY, vZ);
        }

        /* Apply source and save seismogram */
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD3Dacoustic<float>;
template class KITGPI::ForwardSolver::FD3Dacoustic<double>;
