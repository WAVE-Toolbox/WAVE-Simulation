#include "ForwardSolver2Dviscoemem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    return (this->estimateBoundaryMemory(config, dist, modelCoordinates, DampingBoundary, ConvPML));
}

/*! \brief Initialitation of the ForwardSolver
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
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::initForwardSolver(Configuration::Configuration const &config, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Wavefields::WavefieldsEM<ValueType> &wavefield, Modelparameter::ModelparameterEM<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Dviscoemem");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefEX().getDistributionPtr() == model.getMagneticPermeability().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution */
    auto dist = wavefield.getRefEX().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* allocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    
    DT_temp = DT;
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::prepareForModelling(Modelparameter::ModelparameterEM<ValueType> const &model, ValueType DT)
{
    auto const &electricConductivityOpticalAverageX = model.getElectricConductivityOpticalAverageX();
    auto const &electricConductivityOpticalAverageY = model.getElectricConductivityOpticalAverageY();
    auto const &dielectricPermittivityOpticalAverageX = model.getDielectricPermittivityOpticalAverageX();
    auto const &dielectricPermittivityOpticalAverageY = model.getDielectricPermittivityOpticalAverageY();
    auto const &dielectricPermittivityAverageX = model.getDielectricPermittivityAverageX();
    auto const &dielectricPermittivityAverageY = model.getDielectricPermittivityAverageY();
    auto const &tauDielectricPermittivityAverageX = model.getTauDielectricPermittivityAverageX();
    auto const &tauDielectricPermittivityAverageY = model.getTauDielectricPermittivityAverageY();
    ValueType tauElectricDisplacement = model.getTauElectricDisplacement();
    scai::IndexType numRelaxationMechanisms = model.getNumRelaxationMechanisms();

    CaAverageX = this->getAveragedCa(dielectricPermittivityOpticalAverageX, electricConductivityOpticalAverageX, DT);
    CaAverageY = this->getAveragedCa(dielectricPermittivityOpticalAverageY, electricConductivityOpticalAverageY, DT);
    
    CbAverageX = this->getAveragedCb(dielectricPermittivityOpticalAverageX, electricConductivityOpticalAverageX, DT);
    CbAverageY = this->getAveragedCb(dielectricPermittivityOpticalAverageY, electricConductivityOpticalAverageY, DT);
    
    Cc = this->getCc(tauElectricDisplacement, DT);
    
    CdAverageX = this->getAveragedCd(dielectricPermittivityAverageX, tauDielectricPermittivityAverageX, numRelaxationMechanisms, tauElectricDisplacement, DT);
    CdAverageY = this->getAveragedCd(dielectricPermittivityAverageY, tauDielectricPermittivityAverageY, numRelaxationMechanisms, tauElectricDisplacement, DT);
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
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief Running the 2-D visco-emem foward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timesample 
 *
\begin{equation}
\begin{align*}
  &\vec{H}_{z}^{n+\frac{1}{2}} = \vec{H}_{z}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxy}^{-1}\right)  \left(- \underline{D}_{\,x,f}^q \vec{E}_y^{n} + \underline{D}_{\,y,f}^q \vec{E}_x^{n} \right)\\  
  &\vec{E}_{x}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ax}\right) \vec{E}_{x}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{bx}\right) \left[ \frac{1}{\Delta h} \underline{D}_{\,y,b}^q \vec{H}_z^{n+\frac{1}{2}} + \sum_{l=1}^L \vec{r}_{lx}^{n+\frac{1}{2}} \right]\\    
  &\vec{E}_{y}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ay}\right) \vec{E}_{y}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{by}\right) \left[ - \frac{1}{\Delta h} \underline{D}_{\,x,b}^q \vec{H}_z^{n+\frac{1}{2}} + \sum_{l=1}^L \vec{r}_{ly}^{n+\frac{1}{2}} \right]\\   
  &\vec{r}_{lx}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{lx}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l x}\right) \vec{E}_{x}^{n} \\
  &\vec{r}_{ly}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{ly}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l y}\right) \vec{E}_{y}^{n}
\end{align*}
\end{equation}
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Dviscoemem");

    /* Get references to required modelparameter */
    auto const &inverseMagneticPermeabilityAverageXY = model.getInverseMagneticPermeabilityAverageXY();
            
    /* Get references to required wavefields */
    auto &hZ = wavefield.getRefHZ();
    auto &eY = wavefield.getRefEY();
    auto &eX = wavefield.getRefEX();
    auto &rX = wavefield.getRefRX();
    auto &rY = wavefield.getRefRY();

    /* Get references to required derivatives matrixes */
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
    /*  rX rY */
    /* ----------------*/
    // numRelaxationMechanisms = 1
    update = Cc * rX;
    update_temp = CdAverageX * eX;
    rX = update_temp + update;
    
    update = Cc * rY;
    update_temp = CdAverageY * eY;
    rY = update_temp + update;
    
    /* ----------------*/
    /*  ex */
    /* ----------------*/
    update = DybStaggeredX * hZ;
    if (useConvPML) {
        ConvPML.apply_hzy(update);
    }
    update_temp = -DT_temp * rX;
    update += update_temp;
    update *= CbAverageX;    
    update_temp = CaAverageX * eX;
    eX = update_temp + update;
    
    /* ----------------*/
    /*  ey */
    /* ----------------*/
    update_temp = Dxb * hZ;
    if (useConvPML) {
        ConvPML.apply_hzx(update_temp);
    }
    update = -update_temp;
    update_temp = -DT_temp * rY;
    update += update_temp;
    update *= CbAverageY;   
    update_temp = CaAverageY * eY;
    eY = update_temp + update;

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(eY, eX, hZ);
    }

    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscoemem<ValueType>::runAdjoint(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{

}

template class KITGPI::ForwardSolver::FD2Dviscoemem<float>;
template class KITGPI::ForwardSolver::FD2Dviscoemem<double>;
