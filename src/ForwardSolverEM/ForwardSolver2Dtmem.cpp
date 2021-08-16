#include "ForwardSolver2Dtmem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dtmem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
void KITGPI::ForwardSolver::FD2Dtmem<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Dtmem");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefHX().getDistributionPtr() == model.getMagneticPermeability().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefHX().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* allocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    
    DT_temp = DT;
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
void KITGPI::ForwardSolver::FD2Dtmem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dtmem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Dtmem<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    auto const &electricConductivity = model.getElectricConductivity();
    auto const &dielectricPermittivity = model.getDielectricPermittivity();
    
    CaAverageZ = this->getAveragedCa(dielectricPermittivity, electricConductivity, DT);
    CbAverageZ = this->getAveragedCb(dielectricPermittivity, electricConductivity, DT);
}

/*! \brief Running the 2-D tmem foward solver
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
\begin{equation}
\begin{align*}
  &\vec{H}_{x}^{n+\frac{1}{2}} = \vec{H}_{x}^{n-\frac{1}{2}} - \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{myz}^{-1}\right) \underline{D}_{\,y,f}^q \vec{E}_{z}^{n}\\
  &\vec{H}_{y}^{n+\frac{1}{2}} = \vec{H}_{y}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxz}^{-1}\right) \underline{D}_{\,x,f}^q \vec{E}_z^{n} \\
  &\vec{E}_{z}^{n+1} =  \mathrm{diag} \left(\vec{C}_{az}\right) \vec{E}_{z}^{n} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{bz}\right) \left(\underline{D}_{\,x,b}^q \vec{H}_y^{n+\frac{1}{2}} - \underline{D}_{\,y,b}^q \vec{H}_x^{n+\frac{1}{2}} \right)
\end{align*}
\end{equation}
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dtmem<ValueType>::run(Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives)
{
    SCAI_REGION("ForwardSolver.timestep2Dtmem");

    /* Get references to required modelparameter */
    auto const &inverseMagneticPermeabilityAverageYZ = model.getInverseMagneticPermeabilityAverageYZ();
    auto const &inverseMagneticPermeabilityAverageXZ = model.getInverseMagneticPermeabilityAverageXZ();
   
    auto &hX = wavefield.getRefHX();
    auto &hY = wavefield.getRefHY();    
    auto &eZ = wavefield.getRefEZ();
    
    /* Get references to required derivatives matrixes */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();
    auto const &Dyf = derivatives.getDyf();
    auto const &Dyb = derivatives.getDyb();

    /* ----------------*/
    /* hx */
    /* ----------------*/
    update = Dyf * eZ;
    if (useConvPML) {
        ConvPML.apply_ezy(update);
    }
    update *= inverseMagneticPermeabilityAverageYZ;
    hX -= update;

    /* ----------------*/
    /* hy */
    /* ----------------*/
    update_temp = Dxf * eZ;
    if (useConvPML) {
        ConvPML.apply_ezx(update_temp);
    }
    update = -update_temp;
    update *= inverseMagneticPermeabilityAverageXZ;
    hY -= update;
    
    /* ----------------*/
    /*  ez */
    /* ----------------*/
    update = Dxb * hY;
    update_temp = Dyb * hX;
    if (useConvPML) {
        ConvPML.apply_hyx(update);
        ConvPML.apply_hxy(update_temp);
    }
    update -= update_temp;
    update *= CbAverageZ;   
    update_temp = CaAverageZ * eZ;
    eZ = update_temp + update;    

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(eZ, hX, hY);
    }
}

template class KITGPI::ForwardSolver::FD2Dtmem<float>;
template class KITGPI::ForwardSolver::FD2Dtmem<double>;
