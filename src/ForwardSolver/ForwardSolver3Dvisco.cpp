#include "ForwardSolver3Dvisco.hpp"
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
void KITGPI::ForwardSolver::FD3Dvisco<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
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

/*! \brief Running the 3-D visco-elastic foward solver
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
void KITGPI::ForwardSolver::FD3Dvisco<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType tStart, IndexType tEnd, ValueType DT)
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
    lama::Vector const &tauSAverageXY = model.getTauSAverageXY();
    lama::Vector const &tauSAverageXZ = model.getTauSAverageXZ();
    lama::Vector const &tauSAverageYZ = model.getTauSAverageYZ();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getRefVX();
    lama::Vector &vY = wavefield.getRefVY();
    lama::Vector &vZ = wavefield.getRefVZ();

    lama::Vector &Sxx = wavefield.getRefSxx();
    lama::Vector &Syy = wavefield.getRefSyy();
    lama::Vector &Szz = wavefield.getRefSzz();
    lama::Vector &Syz = wavefield.getRefSyz();
    lama::Vector &Sxz = wavefield.getRefSxz();
    lama::Vector &Sxy = wavefield.getRefSxy();

    lama::Vector &Rxx = wavefield.getRefRxx();
    lama::Vector &Ryy = wavefield.getRefRyy();
    lama::Vector &Rzz = wavefield.getRefRzz();
    lama::Vector &Ryz = wavefield.getRefRyz();
    lama::Vector &Rxz = wavefield.getRefRxz();
    lama::Vector &Rxy = wavefield.getRefRxy();

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

    lama::DenseVector<ValueType> update2(vX.getDistributionPtr());

    common::unique_ptr<lama::Vector> vxxPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vyyPtr(vX.newVector());
    common::unique_ptr<lama::Vector> vzzPtr(vX.newVector());

    lama::Vector &vxx = *vxxPtr;
    lama::Vector &vyy = *vyyPtr;
    lama::Vector &vzz = *vzzPtr;

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

        update *= inverseDensityAverageX;
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

        update *= inverseDensityAverageY;
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

        update *= inverseDensityAverageZ;
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
        update *= pWaveModulus;

        update2 = inverseRelaxationTime * update;
        update2 *= tauP;

        Sxx += DThalf * Rxx;
        Rxx *= viscoCoeff1;
        Rxx -= update2;

        Syy += DThalf * Ryy;
        Ryy *= viscoCoeff1;
        Ryy -= update2;

        Szz += DThalf * Rzz;
        Rzz *= viscoCoeff1;
        Rzz -= update2;

        update *= onePlusLtauP;
        Sxx += update;
        Syy += update;
        Szz += update;

        /* Update Sxx and Rxx */
        update = vyy + vzz;
        update *= sWaveModulus;
        update *= 2.0;

        update2 = inverseRelaxationTime * update;
        update2 *= tauS;
        Rxx += update2;
        update *= onePlusLtauS;
        Sxx -= update;

        Rxx *= viscoCoeff2;
        Sxx += DThalf * Rxx;

        /* Update Syy and Ryy */
        update = vxx + vzz;
        update *= sWaveModulus;
        update *= 2.0;

        update2 = inverseRelaxationTime * update;
        update2 *= tauS;
        Ryy += update2;
        update *= onePlusLtauS;
        Syy -= update;

        Ryy *= viscoCoeff2;
        Syy += DThalf * Ryy;

        /* Update Szz and Szz */
        update = vxx + vyy;
        update *= sWaveModulus;
        update *= 2.0;

        update2 = inverseRelaxationTime * update;
        update2 *= tauS;
        Rzz += update2;
        update *= onePlusLtauS;
        Szz -= update;

        Rzz *= viscoCoeff2;
        Szz += DThalf * Rzz;

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

        update *= sWaveModulusAverageXY;

        update2 = inverseRelaxationTime * update;
        update2 *= tauSAverageXY;
        Rxy -= update2;
        update *= onePlusLtauS;
        Sxy += update;

        Rxy *= viscoCoeff2;
        Sxy += DThalf * Rxy;

        /* Update Sxz and Rxz */
        Sxz += DThalf * Rxz;
        Rxz *= viscoCoeff1;

        update = Dzf * vX;
        if (useConvPML) {
            ConvPML.apply_vxz(update);
        }

        update_temp = Dxf * vZ;
        if (useConvPML) {
            ConvPML.apply_vzx(update_temp);
        }
        update += update_temp;

        update *= sWaveModulusAverageXZ;

        update2 = inverseRelaxationTime * update;
        update2 *= tauSAverageXZ;
        Rxz -= update2;
        update *= onePlusLtauS;
        Sxz += update;

        Rxz *= viscoCoeff2;
        Sxz += DThalf * Rxz;

        /* Update Syz and Syz */
        Syz += DThalf * Ryz;
        Ryz *= viscoCoeff1;

        update = Dzf * vY;
        if (useConvPML) {
            ConvPML.apply_vyz(update);
        }

        update_temp = DyfPressure * vZ;
        if (useConvPML) {
            ConvPML.apply_vzy(update_temp);
        }
        update += update_temp;
        update *= sWaveModulusAverageYZ;

        update2 = inverseRelaxationTime * update;
        update2 *= tauSAverageYZ;
        Ryz -= update2;
        update *= onePlusLtauS;
        Syz += update;

        Ryz *= viscoCoeff2;
        Syz += DThalf * Ryz;

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            update = vxx + vzz;
            FreeSurface.apply(update, update2, Sxx, Syy, Szz, Rxx, Ryy, Rzz);
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

template class KITGPI::ForwardSolver::FD3Dvisco<float>;
template class KITGPI::ForwardSolver::FD3Dvisco<double>;
