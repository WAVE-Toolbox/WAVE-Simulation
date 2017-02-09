#include "ForwardSolver3Dacoustic.hpp"
using namespace scai;

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

    /* Prepare Free Surface */
    if (config.get<IndexType>("FreeSurface")) {
        useFreeSurface = true;
        FreeSurface.init(dist, derivatives, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<ValueType>("DH"));
    }

    /* Prepare Damping Boundary */
    if (config.get<IndexType>("DampingBoundary") == 1) {
        if (config.get<IndexType>("DampingBoundaryType") == 1) {
            useDampingBoundary = true;
            DampingBoundary.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"), useFreeSurface);
        }

        if (config.get<IndexType>("DampingBoundaryType") == 2) {
            useConvPML = true;
            ConvPML.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<IndexType>("DH"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("NPower"), config.get<ValueType>("KMaxCPML"), config.get<ValueType>("CenterFrequencyCPML"), config.get<ValueType>("VMaxCPML"), useFreeSurface);
        }
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
 \param NT Total number of time steps
 \param DT Temporal Sampling intervall in seconds
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::run(Acquisition::Receivers<ValueType> &receiver, Acquisition::Sources<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType NT, ValueType /*DT*/)
{

    SCAI_REGION("timestep")

    SCAI_ASSERT_ERROR(NT > 0, " Number of time steps has to be greater than zero. ");

    /* Get references to required modelparameter */
    lama::Vector const &inverseDensity = model.getInverseDensity();
    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector const &inverseDensityAverageZ = model.getInverseDensityAverageZ();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getVX();
    lama::Vector &vY = wavefield.getVY();
    lama::Vector &vZ = wavefield.getVZ();
    lama::Vector &p = wavefield.getP();

    /* Get references to required derivatives matrixes */
    lama::Matrix const &Dxf = derivatives.getDxf();
    lama::Matrix const &Dzf = derivatives.getDzf();
    lama::Matrix const &Dxb = derivatives.getDxb();
    lama::Matrix const &Dzb = derivatives.getDzb();
    lama::Matrix const &Dyb = derivatives.getDyb();
    lama::Matrix const &Dyf = derivatives.getDyfVelocity();

    SourceReceiverImpl::FDTD3Dacoustic<ValueType> SourceReceiver(sources, receiver, wavefield);

    common::unique_ptr<lama::Vector> updatePtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update = *updatePtr;                          // get Reference of VectorPointer

    common::unique_ptr<lama::Vector> update_tempPtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update_temp = *update_tempPtr;                     // get Reference of VectorPointer

    dmemo::CommunicatorPtr comm = inverseDensity.getDistributionPtr()->getCommunicatorPtr();

    /* --------------------------------------- */
    /* Start runtime critical part             */
    /* --------------------------------------- */

    HOST_PRINT(comm, "Start time stepping\n");
    ValueType start_t = common::Walltime::get();
    for (IndexType t = 0; t < NT; t++) {

        if (t % 100 == 0 && t != 0) {
            HOST_PRINT(comm, "Calculating time step " << t << " from " << NT << "\n");
        }

        /* update velocity */
        update = Dxf * p;
        if (useConvPML) {
            ConvPML.apply_p_x(update);
        }
        vX += update.scale(inverseDensityAverageX);

        update = Dyf * p;
        if (useConvPML) {
            ConvPML.apply_p_y(update);
        }
        vY += update.scale(inverseDensityAverageY);

        update = Dzf * p;
        if (useConvPML) {
            ConvPML.apply_p_z(update);
        }
        vZ += update.scale(inverseDensityAverageZ);

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

        p += update.scale(pWaveModulus);

        /* Apply free surface to pressure update */
        if (useFreeSurface) {
            FreeSurface.apply(p);
        }

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(p, vX, vY, vZ);
        }

        /* Apply source and save seismogram */
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);
    }
    ValueType end_t = common::Walltime::get();
    HOST_PRINT(comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}

template class KITGPI::ForwardSolver::FD3Dacoustic<float>;
template class KITGPI::ForwardSolver::FD3Dacoustic<double>;
