#include "ForwardSolver2Demem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Demem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Demem<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Demem");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefHZ().getDistributionPtr() == model.getMagneticPermeability().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefHZ().getDistributionPtr();
    
    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* allocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
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
void KITGPI::ForwardSolver::FD2Demem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Demem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Demem<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    auto const &electricConductivityAverageX = model.getElectricConductivityAverageX();
    auto const &electricConductivityAverageY = model.getElectricConductivityAverageY();
    auto const &dielectricPermittivityAverageX = model.getDielectricPermittivityAverageX();
    auto const &dielectricPermittivityAverageY = model.getDielectricPermittivityAverageY();

    CaAverageX = this->getAveragedCa(dielectricPermittivityAverageX, electricConductivityAverageX, DT);
    CaAverageY = this->getAveragedCa(dielectricPermittivityAverageY, electricConductivityAverageY, DT);
    
    CbAverageX = this->getAveragedCb(dielectricPermittivityAverageX, electricConductivityAverageX, DT);
    CbAverageY = this->getAveragedCb(dielectricPermittivityAverageY, electricConductivityAverageY, DT);
}

/*! \brief Running the 2-Dememforward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timestep
 *
\begin{equation}
\begin{align*}
  &\vec{H}_{z}^{n+\frac{1}{2}} = \vec{H}_{z}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxy}^{-1}\right)  \left(- \underline{D}_{\,x,f}^q \vec{E}_y^{n} + \underline{D}_{\,y,f}^q \vec{E}_x^{n} \right)\\
  &\vec{E}_{x}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ax}\right) \vec{E}_{x}^{n} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{bx}\right) \underline{D}_{\,y,b}^q \vec{H}_z^{n+\frac{1}{2}} \\
  &\vec{E}_{y}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ay}\right) \vec{E}_{y}^{n} - \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{by}\right) \underline{D}_{\,x,b}^q \vec{H}_z^{n+\frac{1}{2}} 
\end{align*}
\end{equation}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Demem<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Demem");
    
    /* Get references to required modelparameter */
    auto const &inverseMagneticPermeabilityAverageXY = model.getInverseMagneticPermeabilityAverageXY();

    /* Get references to required wavefields */
    auto &hZ = wavefield.getRefHZ();
    auto &eY = wavefield.getRefEY();
    auto &eX = wavefield.getRefEX();

    /* Get references to required derivatives matrices */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();   
    auto const &Dyf = derivatives.getDyf();
    
    lama::Matrix<ValueType> const &DybStaggeredX = derivatives.getDybStaggeredX();  

    SourceReceiverImpl::FDTD2Demem<ValueType> SourceReceiver(sources, receiver, wavefield);
    
    /* ----------------*/
    /* hz */
    /* ----------------*/
    update = Dxf * eY;
    update_temp = Dyf * eX;
    if (useConvPML) {
        ConvPML.apply_eyx(update);
        ConvPML.apply_exy(update_temp);
    }
    update -= update_temp;
    update *= inverseMagneticPermeabilityAverageXY;
    hZ -= update;

    /* ----------------*/
    /*  ex */
    /* ----------------*/
    update = DybStaggeredX * hZ;
    if (useConvPML) {
        ConvPML.apply_hzy(update);
    }
    update *= CbAverageX;    
    update_temp = CaAverageX * eX;
    eX = update_temp + update;
    
    /* ----------------*/
    /*  ey */
    /* ----------------*/
    update = Dxb * hZ;
    if (useConvPML) {
        ConvPML.apply_hzx(update);
    }
    update *= CbAverageY;   
    update_temp = CaAverageY * eY;
    eY = update_temp - update;

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(eY, eX, hZ);
    }
    
    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Demem<float>;
template class KITGPI::ForwardSolver::FD2Demem<double>;
