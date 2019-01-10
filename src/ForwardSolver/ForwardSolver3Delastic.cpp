#include "ForwardSolver3Delastic.hpp"
using namespace scai;

/*! \brief Initialitation of the ForwardSolver
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param wavefield Wavefields for the modelling
 \param model model class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType /*DT*/)
{
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefVX().getDistributionPtr();
    ;

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* aalocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    vxx.allocate(dist);
    vyy.allocate(dist);
    vzz.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    vxx.setContextPtr(ctx);
    vyy.setContextPtr(ctx);
    vzz.setContextPtr(ctx);
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/)
{
    if (useFreeSurface == 1) {
        FreeSurface.setModelparameter(model);
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
 *
 * The update equations for velocity, \f$v_i\f$, and stress, \f$\sigma_{ij}\f$, are implemented as follows where \f$M\f$ is the P-wave modulus, \f$\mu\f$ the S-wave modulus and \f$\rho_{inv}\f$ the inverse density. Note that the equations for normal stresses are rearranged compared to the user manual and the scaling with the temporal and spatial discretization is included in the derivative matrices. The velocity update is executed first followed by the stress update and finally the source term is added. If a free surface is chosen, the derivative matrices will be adapted to satisfy the free surface condition.
 *
 \f{eqnarray*}{
	\vec{v}_x &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,f}^q \vec{\sigma}_{xx} + \underline{D}_{\,y,b}^q \vec{\sigma}_{xy} + \underline{D}_{\,z,b}^q \vec{\sigma}_{xz} \right)\\
	\vec{v}_y &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{\sigma}_{yx} + \underline{D}_{\,y,f}^q \vec{\sigma}_{yy} + \underline{D}_{\,z,b}^q \vec{\sigma}_{yz} \right)\\
	\vec{v}_z &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{\sigma}_{zx} + \underline{D}_{\,y,b}^q \vec{\sigma}_{zy} + \underline{D}_{\,z,f}^q \vec{\sigma}_{zz} \right)\\
 \f}
 \f{eqnarray*}{
	\vec{\sigma}_{xx} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y + \underline{D}_{\,z,b}^q \vec{v}_z \right) - 2~ \frac{\Delta t}{\Delta h} ~\mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,y,b}^q \vec{v}_y + \underline{D}_{\,z,b}^q \vec{v}_z \right)\\
	\vec{\sigma}_{yy} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y + \underline{D}_{\,z,b}^q \vec{v}_z \right) - 2~ \frac{\Delta t}{\Delta h} ~\mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x + \underline{D}_{\,z,b}^q \vec{v}_z \right)\\
	\vec{\sigma}_{zz} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{M}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y + \underline{D}_{\,z,b}^q \vec{v}_z \right) - 2~ \frac{\Delta t}{\Delta h} ~\mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{v}_x +\underline{D}_{\,y,b}^q \vec{v}_y \right)
 \f}
 \f{eqnarray*}{
	\vec{\sigma}_{xy} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,y,f}^q \vec{v}_x + \underline{D}_{\,x,f}^q \vec{v}_y \right)\\
	\vec{\sigma}_{xz} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,z,f}^q \vec{v}_x + \underline{D}_{\,x,f}^q \vec{v}_z \right)\\
	\vec{\sigma}_{yz} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,z,f}^q \vec{v}_y + \underline{D}_{\,y,f}^q \vec{v}_z \right)
 \f}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType t)
{

    SCAI_REGION("timestep");

    /* Get references to required modelparameter */
    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();
    lama::Vector<ValueType> const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector<ValueType> const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector<ValueType> const &inverseDensityAverageZ = model.getInverseDensityAverageZ();
    lama::Vector<ValueType> const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();
    lama::Vector<ValueType> const &sWaveModulusAverageXZ = model.getSWaveModulusAverageXZ();
    lama::Vector<ValueType> const &sWaveModulusAverageYZ = model.getSWaveModulusAverageYZ();

    /* Get references to required wavefields */
    lama::Vector<ValueType> &vX = wavefield.getRefVX();
    lama::Vector<ValueType> &vY = wavefield.getRefVY();
    lama::Vector<ValueType> &vZ = wavefield.getRefVZ();

    lama::Vector<ValueType> &Sxx = wavefield.getRefSxx();
    lama::Vector<ValueType> &Syy = wavefield.getRefSyy();
    lama::Vector<ValueType> &Szz = wavefield.getRefSzz();

    lama::Vector<ValueType> &Syz = wavefield.getRefSyz();
    lama::Vector<ValueType> &Sxz = wavefield.getRefSxz();
    lama::Vector<ValueType> &Sxy = wavefield.getRefSxy();

    /* Get references to required derivatives matrixes */
    lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    lama::Matrix<ValueType> const &Dzf = derivatives.getDzf();
    lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    lama::Matrix<ValueType> const &Dzb = derivatives.getDzb();

    lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    lama::Matrix<ValueType> const &DybFreeSurface = derivatives.getDybFreeSurface();
    lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
    lama::Matrix<ValueType> const &DyfFreeSurface = derivatives.getDyfFreeSurface();

    SourceReceiverImpl::FDTD3Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);

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
        update_temp = DybFreeSurface * Sxy;
    } else {
        update_temp = Dyb * Sxy;
    }

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

    update_temp = Dzb * Syz;
    if (useConvPML) {
        ConvPML.apply_syz_z(update_temp);
    }
    update += update_temp;

    update *= inverseDensityAverageY;
    vY += update;

    /* -------- */
    /*    vz    */
    /* -------- */
    update = Dxb * Sxz;
    if (useConvPML) {
        ConvPML.apply_sxz_x(update);
    }

    if (useFreeSurface == 1) {
        /* Apply image method */
        update_temp = DybFreeSurface * Syz;
    } else {
        update_temp = Dyb * Syz;
    }

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

    /* -------------------- */
    /* update normal stress */
    /* -------------------- */
    vxx = Dxb * vX;
    vyy = Dyb * vY;
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

    Sxx += update;
    Syy += update;
    Szz += update;

    update = vyy + vzz;
    update *= sWaveModulus;
    Sxx -= 2.0 * update;
    update = vxx + vzz;
    update *= sWaveModulus;
    Syy -= 2.0 * update;
    update = vxx + vyy;
    update *= sWaveModulus;
    Szz -= 2.0 * update;

    /* ------------------- */
    /* update shear stress */
    /* ------------------- */
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

    update *= sWaveModulusAverageXZ;
    Sxz += update;

    update = Dzf * vY;
    if (useConvPML) {
        ConvPML.apply_vyz(update);
    }
    update_temp = Dyf * vZ;
    if (useConvPML) {
        ConvPML.apply_vzy(update_temp);
    }
    update += update_temp;
    update *= sWaveModulusAverageYZ;
    Syz += update;

    /* Apply free surface to stress update */
    if (useFreeSurface == 1) {
        update = vxx + vzz;
        FreeSurface.exchangeHorizontalUpdate(update, vyy, Sxx, Szz);
    }

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(Sxx, Syy, Szz, Sxy, Sxz, Syz, vX, vY, vZ);
    }

    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD3Delastic<float>;
template class KITGPI::ForwardSolver::FD3Delastic<double>;
