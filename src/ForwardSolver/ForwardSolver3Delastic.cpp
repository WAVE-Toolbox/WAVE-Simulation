#include "ForwardSolver3Delastic.hpp"
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
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

/*! \brief Running the 3-D elastic foward solver
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType tStart, IndexType tEnd, ValueType /*DT*/)
{

    SCAI_REGION("timestep")

    SCAI_ASSERT_ERROR((tEnd - tStart) >= 1, " Number of time steps has to be greater than zero. ");

    /* Get references to required modelparameter */
    lama::Vector const &inverseDensity = model.getInverseDensity();
    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &sWaveModulus = model.getSWaveModulus();
    lama::Vector const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector const &inverseDensityAverageZ = model.getInverseDensityAverageZ();
    lama::Vector const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();
    lama::Vector const &sWaveModulusAverageXZ = model.getSWaveModulusAverageXZ();
    lama::Vector const &sWaveModulusAverageYZ = model.getSWaveModulusAverageYZ();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getVX();
    lama::Vector &vY = wavefield.getVY();
    lama::Vector &vZ = wavefield.getVZ();

    lama::Vector &Sxx = wavefield.getSxx();
    lama::Vector &Syy = wavefield.getSyy();
    lama::Vector &Szz = wavefield.getSzz();

    lama::Vector &Syz = wavefield.getSyz();
    lama::Vector &Sxz = wavefield.getSxz();
    lama::Vector &Sxy = wavefield.getSxy();

    /* Get references to required derivatives matrixes */
    lama::Matrix const &Dxf = derivatives.getDxf();
    lama::Matrix const &Dzf = derivatives.getDzf();
    lama::Matrix const &Dxb = derivatives.getDxb();
    lama::Matrix const &Dzb = derivatives.getDzb();

    lama::Matrix const &DybPressure = derivatives.getDybPressure();
    lama::Matrix const &DybVelocity = derivatives.getDybVelocity();
    lama::Matrix const &DyfPressure = derivatives.getDyfPressure();
    lama::Matrix const &DyfVelocity = derivatives.getDyfVelocity();

    SourceReceiverImpl::FDTD3Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);

    common::unique_ptr<lama::Vector> updatePtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update = *updatePtr;                          // get Reference of VectorPointer

    common::unique_ptr<lama::Vector> update_tempPtr(vX.newVector()); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector &update_temp = *update_tempPtr;                     // get Reference of VectorPointer

    common::unique_ptr<lama::Vector> vxxPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vyyPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vzzPtr(vX.newVector());

    lama::Vector &vxx = *vxxPtr;
    lama::Vector &vyy = *vyyPtr;
    lama::Vector &vzz = *vzzPtr;

    if (useFreeSurface) {
        FreeSurface.setModelparameter(model);
    }

    dmemo::CommunicatorPtr comm = inverseDensity.getDistributionPtr()->getCommunicatorPtr();

    /* --------------------------------------- */
    /* Start runtime critical part             */
    /* --------------------------------------- */

    for (IndexType t = tStart; t < tEnd; t++) {

        if (t % 100 == 0 && t != 0) {
            HOST_PRINT(comm, "Calculating time step " << t << "\n");
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

        update_temp = Dzb * Sxz;
        if (useConvPML) {
            ConvPML.apply_sxz_z(update_temp);
        }
        update += update_temp;
        update.scale(inverseDensityAverageX);
        vX += update;

        update = Dxb * Sxy;
        if (useConvPML) {
            ConvPML.apply_sxy_x(update);
        }

        update_temp = DyfVelocity * Syy;
        if (useConvPML) {
            ConvPML.apply_syy_y(update_temp);
        }
        update += update_temp;

        update_temp = Dzb * Syz;
        if (useConvPML) {
            ConvPML.apply_syz_z(update_temp);
        }
        update += update_temp;

        update.scale(inverseDensityAverageY);
        vY += update;

        update = Dxb * Sxz;
        if (useConvPML) {
            ConvPML.apply_sxz_x(update);
        }

        update_temp = DybVelocity * Syz;
        if (useConvPML) {
            ConvPML.apply_syz_y(update_temp);
        }
        update += update_temp;

        update_temp = Dzf * Szz;
        if (useConvPML) {
            ConvPML.apply_szz_z(update_temp);
        }
        update += update_temp;

        update.scale(inverseDensityAverageZ);
        vZ += update;

        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        vzz = Dzb * vZ;
        if (useConvPML) {
            ConvPML.apply_vxx(vxx);
            ConvPML.apply_vyy(vyy);
            ConvPML.apply_vzz(vzz);
        }

        update = vxx;
        update += vyy;
        update += vzz;
        update.scale(pWaveModulus);

        Sxx += update;
        Syy += update;
        Szz += update;

        update = vyy + vzz;
        update.scale(sWaveModulus);
        Sxx -= 2.0 * update;
        update = vxx + vzz;
        update.scale(sWaveModulus);
        Syy -= 2.0 * update;
        update = vxx + vyy;
        update.scale(sWaveModulus);
        Szz -= 2.0 * update;

        update = DyfPressure * vX;
        if (useConvPML) {
            ConvPML.apply_vxy(update);
        }
        update_temp = Dxf * vY;
        if (useConvPML) {
            ConvPML.apply_vyx(update_temp);
        }
        update += update_temp;
        update.scale(sWaveModulusAverageXY);
        Sxy += update;

        update = Dzf * vX;
        if (useConvPML) {
            ConvPML.apply_vxz(update);
        }
        update_temp = Dxf * vZ;
        if (useConvPML) {
            ConvPML.apply_vzx(update_temp);
        }
        update += update_temp;

        update.scale(sWaveModulusAverageXZ);
        Sxz += update;

        update = Dzf * vY;
        if (useConvPML) {
            ConvPML.apply_vyz(update);
        }
        update_temp = DyfPressure * vZ;
        if (useConvPML) {
            ConvPML.apply_vzy(update_temp);
        }
        update += update_temp;
        update.scale(sWaveModulusAverageYZ);
        Syz += update;

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            update = vxx + vzz;
            FreeSurface.apply(update, Sxx, Syy, Szz);
        }

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(Sxx, Syy, Szz, Sxy, Sxz, Syz, vX, vY, vZ);
        }

        /* Apply source and save seismogram */
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);
    }

    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}

template class KITGPI::ForwardSolver::FD3Delastic<float>;
template class KITGPI::ForwardSolver::FD3Delastic<double>;