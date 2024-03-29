#include "ForwardSolver2Dviscoelastic.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    return (this->estimateBoundaryMemory(config, dist, modelCoordinates, DampingBoundary, ConvPML));
}

/*! \brief Initialization of the ForwardSolver
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param wavefield Wavefields for the modelling
 \param model model class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param DT time sampling interval
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Dviscoelastic");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution */
    auto dist = wavefield.getRefVX().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* allocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    vxx.allocate(dist);
    vyy.allocate(dist);
    update2.allocate(dist);
    onePlusLtauP.allocate(dist);
    onePlusLtauS.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    vxx.setContextPtr(ctx);
    vyy.setContextPtr(ctx);
    update2.setContextPtr(ctx);
    onePlusLtauP.setContextPtr(ctx);
    onePlusLtauS.setContextPtr(ctx);

    numRelaxationMechanisms = model.getNumRelaxationMechanisms();
    viscoCoeff1.clear();
    viscoCoeff2.clear();
    relaxationTime.clear();
    inverseRelaxationTime.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * model.getRelaxationFrequency()[l])); // = 1 / ( 2 * Pi * f_relax )
        inverseRelaxationTime.push_back(1.0 / relaxationTime[l]);
        viscoCoeff1.push_back(1.0 - DT / (2.0 * relaxationTime[l]));
        viscoCoeff2.push_back(1.0 / (1.0 + DT / (2.0 * relaxationTime[l])));              // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
    }
    DThalf = DT / 2.0;
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::resetCPML()
{
    if (useConvPML) {
        ConvPML.resetCPML();
    }
}

/*! \brief Preparations before each modelling
 *
 *
 \param model model parameter object

 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    /* Get reference to required model vectors */
    lama::Vector<ValueType> const &tauS = model.getTauS();
    lama::Vector<ValueType> const &tauP = model.getTauP();

    onePlusLtauP = 1.0;
    onePlusLtauP += numRelaxationMechanisms * tauP;

    onePlusLtauS = 1.0;
    onePlusLtauS += numRelaxationMechanisms * tauS;

    if (useFreeSurface == 1) {
        FreeSurface.setModelparameter(model, onePlusLtauP, onePlusLtauS, DT);
    }
}

/*! \brief Initialization of the boundary conditions
 *
 *
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief Running the 2-D viscoelastic forward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timesample 
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoelastic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Dviscoelastic");

    /* Get references to required modelparameter */
    auto const &pWaveModulus = model.getPWaveModulus();
    auto const &sWaveModulus = model.getSWaveModulus();
    auto const &inverseDensityAverageX = model.getInverseDensityAverageX();
    auto const &inverseDensityAverageY = model.getInverseDensityAverageY();
    auto const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();
    auto const &tauSAverageXY = model.getTauSAverageXY();

    /* Get references to required wavefields */
    auto &vX = wavefield.getRefVX();
    auto &vY = wavefield.getRefVY();

    auto &Sxx = wavefield.getRefSxx();
    auto &Syy = wavefield.getRefSyy();
    auto &Sxy = wavefield.getRefSxy();

    auto &Rxx = wavefield.getRefRxx();
    auto &Ryy = wavefield.getRefRyy();
    auto &Rxy = wavefield.getRefRxy();

    /* Get references to required derivatives matrices */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();

    auto const &Dyb = derivatives.getDyb();
    auto const &DybFreeSurface = derivatives.getDybFreeSurface();
    auto const &Dyf = derivatives.getDyf();
    auto const &DyfFreeSurface = derivatives.getDyfFreeSurface(); 

    /* Get reference to required model vectors */
    auto const &tauS = model.getTauS();
    auto const &tauP = model.getTauP();
    
    SourceReceiverImpl::FDTD2Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);
    
    /* ----------------*/
    /* update velocity */
    /* ----------------*/
    update = Dxf * Sxx;
    if (useConvPML) {
        ConvPML.apply_sxx_x(update);
    }

    if (useFreeSurface == 1) {
        /* Apply image method */
        update_temp = DybFreeSurface * Sxy;
    } else {
        update_temp = Dyb * Sxy;
    }
    if (useConvPML) {
        ConvPML.apply_sxy_y(update_temp);
    }
    update += update_temp;
    update *= inverseDensityAverageX;
    vX += update;

    update = Dxb * Sxy;
    if (useConvPML) {
        ConvPML.apply_sxy_x(update);
    }

    if (useFreeSurface == 1) {
        /* Apply image method */
        update_temp = DyfFreeSurface * Syy;
    } else {
        update_temp = Dyf * Syy;
    }
    if (useConvPML) {
        ConvPML.apply_syy_y(update_temp);
    }
    update += update_temp;

    update *= inverseDensityAverageY;
    vY += update;

    /* ----------------*/
    /* pressure update */
    /* ----------------*/

    vxx = Dxb * vX;
    vyy = Dyb * vY;
    if (useConvPML) {
        ConvPML.apply_vxx(vxx);
        ConvPML.apply_vyy(vyy);
    }

    update = vxx;
    update += vyy;
    update *= pWaveModulus;

    for (int l=0; l<numRelaxationMechanisms; l++) {
        update2 = inverseRelaxationTime[l] * update;
        update2 *= tauP;

        Sxx += DThalf * Rxx[l];
        Rxx[l] *= viscoCoeff1[l];
        Rxx[l] -= update2;

        Syy += DThalf * Ryy[l];
        Ryy[l] *= viscoCoeff1[l];
        Ryy[l] -= update2;
    }
    update *= onePlusLtauP;
    Sxx += update;
    Syy += update;

    /* Update Sxx and Rxx */
    update = vyy;
    update *= sWaveModulus;
    update *= 2.0;

    for (int l=0; l<numRelaxationMechanisms; l++) {
        update2 = inverseRelaxationTime[l] * update;

        update2 *= tauS;
        Rxx[l] += update2;

        Rxx[l] *= viscoCoeff2[l];
        Sxx += DThalf * Rxx[l];
    }
    update *= onePlusLtauS;
    Sxx -= update;

    /* Update Syy and Ryy */
    update = vxx;
    update *= sWaveModulus;
    update *= 2.0;

    for (int l=0; l<numRelaxationMechanisms; l++) {
        update2 = inverseRelaxationTime[l] * update;
        update2 *= tauS;
        Ryy[l] += update2;

        Ryy[l] *= viscoCoeff2[l];
        Syy += DThalf * Ryy[l];
    }
    update *= onePlusLtauS;
    Syy -= update;

    /* Update Sxy and Rxy*/
    update = Dyf * vX;
    if (useConvPML) {
        ConvPML.apply_vxy(update);
    }

    update_temp = Dxf * vY;
    if (useConvPML) {
        ConvPML.apply_vyx(update_temp);
    }
    update += update_temp;

    update *= sWaveModulusAverageXY;
    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Sxy += DThalf * Rxy[l];
        Rxy[l] *= viscoCoeff1[l];

        update2 = inverseRelaxationTime[l] * update;
        update2 *= tauSAverageXY;
        Rxy[l] -= update2;

        Rxy[l] *= viscoCoeff2[l];
        Sxy += DThalf * Rxy[l];
    }
    update *= onePlusLtauS;
    Sxy += update;

    /* Apply free surface to stress update */
    if (useFreeSurface) {
        FreeSurface.exchangeHorizontalUpdate(vxx, vyy, Sxx, Rxx, DThalf);
        FreeSurface.setSurfaceZero(Syy);
        for (int l=0; l<numRelaxationMechanisms; l++) {
            FreeSurface.setSurfaceZero(Ryy[l]);
        }
    }

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(Sxx, Syy, Sxy, vX, vY);
    }
    
    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Dviscoelastic<float>;
template class KITGPI::ForwardSolver::FD2Dviscoelastic<double>;
