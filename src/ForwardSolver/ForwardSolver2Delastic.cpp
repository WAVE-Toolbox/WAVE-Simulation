#include "ForwardSolver2Delastic.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Delastic<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    return (this->estimateBoundaryMemory(config, dist, modelCoordinates, DampingBoundary, ConvPML));
}

/*! \brief Initialitation of the ForwardSolver
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param wavefield Wavefields for the modelling
 \param model
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/)
{
    SCAI_REGION("ForwardSolver.init2Delastic");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefVX().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* aalocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    vxx.allocate(dist);
    vyy.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    vxx.setContextPtr(ctx);
    vyy.setContextPtr(ctx);
}

/*! \brief Initialitation of the boundary conditions (wrapper for prepareBoundaries in forwardsolver.cpp)
 *
 *
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/)
{
    if (useFreeSurface == 1) {
        FreeSurface.setModelparameter(model);
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
 \param tStart Counter start in for loop over time steps
 \param tEnd Counter end  in for loop over time steps
 \param DT Temporal Sampling intervall in seconds
  *
  * The update equations for velocity, \f$v_i\f$, and stress, \f$\sigma_{ij}\f$, are implemented as follows where \f$M\f$ is the P-wave modulus, \f$\mu\f$ the S-wave modulus and \f$\rho_{inv}\f$ the inverse density. Note that the equations for normal stresses are rearranged compared to the user manual and the scaling with the temporal and spatial discretization is included in the derivative matrices. The velocity update is executed first followed by the stress update and finally the source term is added. If a free surface is chosen, the derivative matrices will be adapted to satisfy the free surface condition.
 *
 \f{eqnarray*}{
	\vec{v}_x &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,f}^q \vec{\sigma}_{xx} + \underline{D}_{\,y,b}^q \vec{\sigma}_{xy} \right)\\
	\vec{v}_y &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{\sigma}_{yx} + \underline{D}_{\,y,f}^q \vec{\sigma}_{yy} \right)
 \f}
 \f{eqnarray*}{
	\vec{\sigma}_{xx} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y \right) - 2~ \frac{\Delta t}{\Delta h} ~\mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,y,b}^q \vec{v}_y \right)\\
	\vec{\sigma}_{yy} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y \right) - 2~ \frac{\Delta t}{\Delta h} ~\mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x \right)
 \f} 
 \f{eqnarray*}{
	\vec{\sigma}_{xy} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,y,f}^q \vec{v}_x + \underline{D}_{\,x,f}^q \vec{v}_y \right)
 \f}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Delastic");

    /* Get references to required modelparameter */
    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();
    lama::Vector<ValueType> const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector<ValueType> const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector<ValueType> const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();

    /* Get references to required wavefields */
    lama::Vector<ValueType> &vX = wavefield.getRefVX();
    lama::Vector<ValueType> &vY = wavefield.getRefVY();

    lama::Vector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::Vector<ValueType> &Syy = wavefield.getRefSyy();

    lama::Vector<ValueType> &Sxy = wavefield.getRefSxy();

    /* Get references to required derivatives matrixes */
    lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();

    lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    lama::Matrix<ValueType> const &DybStaggeredX = derivatives.getDybStaggeredX();
    lama::Matrix<ValueType> const &DybStaggeredXFreeSurface = derivatives.getDybStaggeredXFreeSurface();

    //   lama::Matrix<ValueType> const &DybFreeSurface = derivatives.getDybFreeSurface();
    lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
    lama::Matrix<ValueType> const &DyfStaggeredX = derivatives.getDyfStaggeredX();
    lama::Matrix<ValueType> const &DyfFreeSurface = derivatives.getDyfFreeSurface();

    /* Get pointers to required interpolation matrices (optional) */
    lama::Matrix<ValueType> const *DinterpolateFull = derivatives.getInterFull();
    lama::Matrix<ValueType> const *DinterpolateStaggeredX = derivatives.getInterStaggeredX();

    SourceReceiverImpl::FDTD2Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);

    /* ----------------*/
    /* update velocity */
    /* ----------------*/

    /* -------- */
    /*    vx    */
    /* -------- */
    update = Dxf * Sxx;
    if (useConvPML) {
        ConvPML.apply_sxx_x(update);
    }

    if (useFreeSurface == 1) {
        /* Apply image method */
        update_temp = DybStaggeredXFreeSurface * Sxy;
    } else {
        update_temp = DybStaggeredX * Sxy;
    }

    if (useConvPML) {
        ConvPML.apply_sxy_y(update_temp);
    }
    update += update_temp;

    update *= inverseDensityAverageX;
    vX += update;

    if (DinterpolateStaggeredX) {
        /* interpolation for vx ghost points at the variable grid interfaces*/
        update_temp.swap(vX);
        vX = *DinterpolateStaggeredX * update_temp;
    }

    /* -------- */
    /*    vy    */
    /* -------- */
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

    if (DinterpolateFull) {
        /* interpolation for vy ghost pointsa t the variable grid interfaces.
         This interpolation has no effect on the simulation.
         Nethertheless it will be done to avoid abitrary values.
         This is helpful for applications like FWI*/
        update_temp.swap(vY);
        vY = *DinterpolateFull * update_temp;
    }

    /* -------------------- */
    /* update normal stress */
    /* -------------------- */
    vxx = Dxb * vX;
    vyy = Dyb * vY;
    if (useConvPML) {
        ConvPML.apply_vxx(vxx);
        ConvPML.apply_vyy(vyy);
    }

    update = vxx;
    update += vyy;
    update *= pWaveModulus;

    Sxx += update;
    Syy += update;

    update = vyy;
    update *= sWaveModulus;
    Sxx -= 2.0 * update;
    update = vxx;
    update *= sWaveModulus;
    Syy -= 2.0 * update;

    if (DinterpolateFull) {
        // interpolation for Sxx/Sxx ghost points at the variable grid interfaces.
        update_temp.swap(Sxx);
        Sxx = *DinterpolateFull * update_temp;

        update_temp.swap(Syy);
        Syy = *DinterpolateFull * update_temp;
    }

    /* ------------------- */
    /* update shear stress */
    /* ------------------- */
    update = DyfStaggeredX * vX;
    if (useConvPML) {
        ConvPML.apply_vxy(update);
    }

    update_temp = Dxf * vY;
    if (useConvPML) {
        ConvPML.apply_vyx(update_temp);
    }
    update += update_temp;

    update *= sWaveModulusAverageXY;
    Sxy += update;

    if (DinterpolateStaggeredX) {
        /* interpolation for Sxy ghost points at the variable grid interfaces.
         This interpolation has no effect on the simulation.
         Nethertheless it will be done to avoid abitrary values.
         This is helpful for applications like FWI*/
        update_temp.swap(Sxy);
        Sxy = *DinterpolateStaggeredX * update_temp;
    }

    /* Apply free surface to horizontal stress update */
    if (useFreeSurface == 1) {
        FreeSurface.exchangeHorizontalUpdate(vxx, vyy, Sxx);
        FreeSurface.setSurfaceZero(Syy);
    }

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(Sxx, Syy, Sxy, vX, vY);
    }

    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Delastic<float>;
template class KITGPI::ForwardSolver::FD2Delastic<double>;
