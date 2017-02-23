#include "ForwardSolver2Dvisco.hpp"
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
void KITGPI::ForwardSolver::FD2Dvisco<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
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

/*! \brief Running the 2-D visco-elastic foward solver
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
void KITGPI::ForwardSolver::FD2Dvisco<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType NT, ValueType DT)
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
    lama::Vector const &tauSAverageXY = model.getTauSAverageXY();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getVX();
    lama::Vector &vY = wavefield.getVY();

    lama::Vector &Sxx = wavefield.getSxx();
    lama::Vector &Syy = wavefield.getSyy();
    lama::Vector &Sxy = wavefield.getSxy();

    lama::Vector &Rxx = wavefield.getRxx();
    lama::Vector &Ryy = wavefield.getRyy();
    lama::Vector &Rxy = wavefield.getRxy();

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

    lama::DenseVector<ValueType> update2(vX.getDistributionPtr());

    common::unique_ptr<lama::Vector> vxxPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vyyPtr(vX.newVector());

    lama::Vector &vxx = *vxxPtr;
    lama::Vector &vyy = *vyyPtr;

    lama::Vector const &tauS = model.getTauS();
    lama::Vector const &tauP = model.getTauP();

    IndexType numRelaxationMechanisms = model.getNumRelaxationMechanisms();         // = Number of relaxation mechanisms
    ValueType relaxationTime = 1.0 / (2.0 * M_PI * model.getRelaxationFrequency()); // = 1 / ( 2 * Pi * f_relax )
    ValueType inverseRelaxationTime = 1.0 / relaxationTime;                         // = 1 / relaxationTime
    ValueType viscoCoeff1 = (1.0 - DT / (2.0 * relaxationTime));                    // = 1 - DT / ( 2 * tau_Sigma_l )
    ValueType viscoCoeff2 = 1.0 / (1.0 + DT / (2.0 * relaxationTime));              // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
    ValueType DThalf = DT / 2.0;                                                    // = DT / 2.0

    lama::DenseVector<ValueType> onePlusLtauP(vX.getDistributionPtr()); // = ( 1 + L * tauP )
    lama::DenseVector<ValueType> onePlusLtauS(vX.getDistributionPtr()); // = ( 1 + L * tauS )

    onePlusLtauP = 1.0;
    onePlusLtauP += numRelaxationMechanisms * tauP;

    onePlusLtauS = 1.0;
    onePlusLtauS += numRelaxationMechanisms * tauS;

    if (useFreeSurface) {
        FreeSurface.setModelparameter(model, onePlusLtauP, onePlusLtauS);
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

	update.scale(inverseDensityAverageY);
        vY += update;

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

        update2 = inverseRelaxationTime * update;
        update2.scale(tauP);

        Sxx += DThalf * Rxx;
        Rxx *= viscoCoeff1;
        Rxx -= update2;

        Syy += DThalf * Ryy;
        Ryy *= viscoCoeff1;
        Ryy -= update2;

        update.scale(onePlusLtauP);
        Sxx += update;
        Syy += update;

        /* Update Sxx and Rxx */
        vyy.scale(sWaveModulus);
        vyy *= 2.0;

        update2 = inverseRelaxationTime * vyy;
	
	update2.scale(tauS);
        Rxx += update2;
	vyy.scale(onePlusLtauS);
        Sxx -= vyy;

        Rxx *= viscoCoeff2;
        Sxx += DThalf * Rxx;

        /* Update Syy and Ryy */
        vxx.scale(sWaveModulus);
        vxx *= 2.0;

        update2 = inverseRelaxationTime * vxx;
	update2.scale(tauS);
        Ryy += update2;
	vxx.scale(onePlusLtauS);
        Syy -= vxx;

        Ryy *= viscoCoeff2;
        Syy += DThalf * Ryy;

        /* Update Sxy and Rxy*/
        Sxy += DThalf * Rxy;
        Rxy *= viscoCoeff1;

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

        update2 = inverseRelaxationTime * update;
	update2.scale(tauSAverageXY);
        Rxy -= update2;
	update.scale(onePlusLtauS);
        Sxy += update;

        Rxy *= viscoCoeff2;
        Sxy += DThalf * Rxy;

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            FreeSurface.apply(vxx, update2, Sxx, Syy, Rxx, Ryy);
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

template class KITGPI::ForwardSolver::FD2Dvisco<float>;
template class KITGPI::ForwardSolver::FD2Dvisco<double>;
