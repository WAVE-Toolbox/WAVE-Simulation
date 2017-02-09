#include "ForwardSolver2Delastic.hpp"
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
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

/*! \brief Running the 2-D elastic foward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::run(Acquisition::Receivers<ValueType> &receiver, Acquisition::Sources<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType NT, ValueType /*DT*/)
{

    SCAI_REGION("timestep")

    SCAI_ASSERT_ERROR(NT > 0, " Number of time steps has to be greater than zero. ");

    /* Get references to required modelparameter */
    lama::Vector const &inverseDensity = model.getInverseDensity();
    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &sWaveModulus = model.getSWaveModulus();
    lama::Vector const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getVX();
    lama::Vector &vY = wavefield.getVY();

    lama::Vector &Sxx = wavefield.getSxx();
    lama::Vector &Syy = wavefield.getSyy();

    lama::Vector &Sxy = wavefield.getSxy();

    /* Get references to required derivatives matrixes */
    lama::Matrix const &Dxf = derivatives.getDxf();
    lama::Matrix const &Dxb = derivatives.getDxb();

    lama::Matrix const &DybPressure = derivatives.getDybPressure();
    lama::Matrix const &DybVelocity = derivatives.getDybVelocity();
    lama::Matrix const &DyfPressure = derivatives.getDyfPressure();
    lama::Matrix const &DyfVelocity = derivatives.getDyfVelocity();

    SourceReceiverImpl::FDTD2Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);

    common::unique_ptr<lama::Vector> updatePtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update = *updatePtr;                          // get Reference of VectorPointer

    common::unique_ptr<lama::Vector> update_tempPtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update_temp = *update_tempPtr;                     // get Reference of VectorPointer

    common::unique_ptr<lama::Vector> vxxPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vyyPtr(vX.newVector());

    lama::Vector &vxx = *vxxPtr;
    lama::Vector &vyy = *vyyPtr;

    if (useFreeSurface) {
        FreeSurface.setModelparameter(model);
    }

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

        /* ----------------*/
        /* update velocity */
        /* ----------------*/
        update = Dxf * Sxx;
        if (useConvPML) {
            ConvPML.apply_sxx_x(update);
        }

        update_temp = DybVelocity * Sxy;
        if (useConvPML) {
            ConvPML.apply_sxy_y(update_temp);
        }
        update += update_temp;

        vX += update.scale(inverseDensityAverageX);

        update = Dxb * Sxy;
        if (useConvPML) {
            ConvPML.apply_sxy_x(update);
        }
        update_temp = DyfVelocity * Syy;
        if (useConvPML) {
            ConvPML.apply_syy_y(update_temp);
        }
        update += update_temp;

        vY += update.scale(inverseDensityAverageY);

        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        if (useConvPML) {
            ConvPML.apply_vxx(vxx);
            ConvPML.apply_vyy(vyy);
        }

        update = vxx;
        update += vyy;
        update.scale(pWaveModulus);

        Sxx += update;
        Syy += update;

        Sxx -= 2.0 * vyy.scale(sWaveModulus);
        Syy -= 2.0 * vxx.scale(sWaveModulus);

        update = DyfPressure * vX;
        if (useConvPML) {
            ConvPML.apply_vxy(update);
        }

        update_temp = Dxf * vY;
        if (useConvPML) {
            ConvPML.apply_vyx(update_temp);
        }
        update += update_temp;

        Sxy += update.scale(sWaveModulusAverageXY);

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            FreeSurface.apply(vxx, Sxx, Syy);
        }

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(Sxx, Syy, Sxy, vX, vY);
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

template class KITGPI::ForwardSolver::FD2Delastic<float>;
template class KITGPI::ForwardSolver::FD2Delastic<double>;
