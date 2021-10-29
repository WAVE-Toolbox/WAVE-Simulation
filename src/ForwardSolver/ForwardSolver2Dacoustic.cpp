#include "ForwardSolver2Dacoustic.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dacoustic<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
void KITGPI::ForwardSolver::FD2Dacoustic<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/)
{
    SCAI_REGION("ForwardSolver.init2Dacoustic");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefVX().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* Initialisation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
}

/*! \brief Initialitation of the boundary conditions
 *
 *
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dacoustic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dacoustic<ValueType>::resetCPML()
{
    if (useConvPML) {
        ConvPML.resetCPML();
    }
}

/*! \brief Running the 2-D acoustic foward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current time sample
 *
 * The update equations for velocity, \f$v_i\f$, and pressure, \f$p\f$, are implemented as follows where \f$M\f$ is the P-wave modulus and \f$\rho_{inv}\f$ the inverse density. Note that the scaling with the temporal and spatial discretization is included in the derivative matrices. The velocity update is executed first followed by the pressure update and finally the source term is added. If a free surface is chosen, the derivative matrices will be adapted to satisfy the free surface condition.
 *
 \f{eqnarray*}
	\vec{v}_x &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,f}^q \vec{p} \right)\\
	\vec{v}_y &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,y,f}^q \vec{p} \right)
 \f}
 \f{eqnarray*}
	\vec{p} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y \right) 
 \f}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dacoustic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Dacoustic");

    /* Get references to required modelparameter */
    auto const &pWaveModulus = model.getPWaveModulus();
    auto const &inverseDensityAverageX = model.getInverseDensityAverageX();
    auto const &inverseDensityAverageY = model.getInverseDensityAverageY();

    /* Get references to required wavefields */
    auto &vX = wavefield.getRefVX();
    auto &vY = wavefield.getRefVY();
    auto &p = wavefield.getRefP();

    /* Get references to required derivatives matrixes */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();
    auto const &Dyb = derivatives.getDyb();
    auto const &Dyf = derivatives.getDyf();
    auto const &DyfFreeSurface = derivatives.getDyfFreeSurface();

    /* Get pointers to required interpolation matrices (optional) */
    auto const *DinterpolateFull = derivatives.getInterFull();
    auto const *DinterpolateStaggeredX = derivatives.getInterStaggeredX();

    SourceReceiverImpl::FDTD2Dacoustic<ValueType> SourceReceiver(sources, receiver, wavefield);

    /* ----------------*/
    /* update velocity */
    /* ----------------*/    
    update = Dxf * p;
    if (useConvPML) {
        ConvPML.apply_p_x(update);
    }
    update *= inverseDensityAverageX;
    vX += update;

    if (DinterpolateStaggeredX) {
        /* interpolation for vz ghost points at the variable grid interfaces
        This interpolation has no effect on the simulation.
         Nevertheless it will be done to avoid arbitrary values.
         This is helpful for applications like FWI*/
        update_temp.swap(vX);
        vX = *DinterpolateStaggeredX * update_temp;
    }

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

    if (DinterpolateFull) {
        /* interpolation for vz ghost points at the variable grid interfaces
        This interpolation has no effect on the simulation.
         Nevertheless it will be done to avoid arbitrary values.
         This is helpful for applications like FWI*/
        update_temp.swap(vY);
        vY = *DinterpolateFull * update_temp;
    }

    /* --------------- */
    /* update pressure */
    /* --------------- */
    update = Dxb * vX;
    if (useConvPML) {
        ConvPML.apply_vxx(update);
    }

    update_temp = Dyb * vY;
    if (useConvPML) {
        ConvPML.apply_vyy(update_temp);
    }
    update += update_temp;

    update *= pWaveModulus;
    p += update;

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(p, vX, vY);
    }

    if (DinterpolateFull) {
        // interpolation for missing pressure points
        update_temp.swap(p);
        p = *DinterpolateFull * update_temp;
    }

    if (useFreeSurface == 1) {
        FreeSurface.setSurfaceZero(p);
    }
    
    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Dacoustic<float>;
template class KITGPI::ForwardSolver::FD2Dacoustic<double>;
