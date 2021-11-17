#include "ForwardSolver2Dviscotmem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    return (this->estimateBoundaryMemory(config, dist, modelCoordinates, DampingBoundary, ConvPML));
}

/*! \brief Initialization of the ForwardSolver
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
void KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Dviscotmem");
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

/*! \brief Initialization of the boundary conditions (wrapper for prepareBoundaries in forwardsolver.cpp)
 *
 *
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    scai::lama::DenseVector<ValueType> electricConductivityEffectiveOptical;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEffectiveOptical;
    auto const &dielectricPermittivity = model.getDielectricPermittivity();
    auto const &electricConductivity = model.getElectricConductivity();
    auto const &tauDielectricPermittivity = model.getTauDielectricPermittivity();
    auto const &tauElectricConductivity = model.getTauElectricConductivity();
    
    numRelaxationMechanisms = model.getNumRelaxationMechanisms();
    std::vector<ValueType> relaxationTime;          // = 1 / ( 2 * Pi * f_relax )
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * model.getRelaxationFrequency()[l])); // = 1 / ( 2 * Pi * f_relax )
    }
    
    dielectricPermittivityEffectiveOptical = this->getDielectricPermittivityEffectiveOptical(dielectricPermittivity, electricConductivity, tauDielectricPermittivity, tauElectricConductivity, model.getDielectricPermittivityVacuum());
    electricConductivityEffectiveOptical = this->getElectricConductivityEffectiveOptical(dielectricPermittivity, electricConductivity, tauDielectricPermittivity, numRelaxationMechanisms, relaxationTime);

    CaAverageZ = this->getAveragedCa(dielectricPermittivityEffectiveOptical, electricConductivityEffectiveOptical, DT);
    
    CbAverageZ = this->getAveragedCb(dielectricPermittivityEffectiveOptical, electricConductivityEffectiveOptical, DT);
        
    Cc = this->getCc(numRelaxationMechanisms, relaxationTime, DT);
    
    CdAverageZ = this->getAveragedCd(dielectricPermittivity, tauDielectricPermittivity, numRelaxationMechanisms, relaxationTime, DT);
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
  &\vec{H}_{x}^{n+\frac{1}{2}} = \vec{H}_{x}^{n-\frac{1}{2}} - \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{myz}^{-1}\right) \underline{D}_{\,y,f}^q \vec{E}_{z}^{n} \\
  &\vec{H}_{y}^{n+\frac{1}{2}} = \vec{H}_{y}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxz}^{-1}\right) \underline{D}_{\,x,f}^q \vec{E}_z^{n} \\
  &\vec{E}_{z}^{n+1} =  \mathrm{diag} \left(\vec{C}_{az}\right) \vec{E}_{z}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{bz}\right) \left[ \frac{1}{\Delta h} \left(\underline{D}_{\,x,b}^q \vec{H}_y^{n+\frac{1}{2}} - \underline{D}_{\,y,b}^q \vec{H}_x^{n+\frac{1}{2}} \right) + \sum_{l=1}^L \vec{r}_{lz}^{n+\frac{1}{2}} \right]\\ 
  &\vec{r}_{lz}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{lz}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l z}\right) \vec{E}_{z}^{n}
\end{align*}
\end{equation}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscotmem<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Dtmem");

    /* Get references to required modelparameter */
    auto const &inverseMagneticPermeabilityAverageYZ = model.getInverseMagneticPermeabilityAverageYZ();
    auto const &inverseMagneticPermeabilityAverageXZ = model.getInverseMagneticPermeabilityAverageXZ();
        
    /* Get references to required wavefields */
    auto &hX = wavefield.getRefHX();
    auto &hY = wavefield.getRefHY();    
    auto &eZ = wavefield.getRefEZ();  
    auto &rZ = wavefield.getRefRZ();

    /* Get references to required derivatives matrices */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();
    auto const &Dyf = derivatives.getDyf();
    auto const &Dyb = derivatives.getDyb();

    SourceReceiverImpl::FDTD2Dtmem<ValueType> SourceReceiver(sources, receiver, wavefield);
    
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
    /*  rZ             */
    /* ----------------*/    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        update = Cc[l] * rZ[l];
        update_temp = CdAverageZ[l] * eZ;
        rZ[l] = update_temp + update;
    }
        
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
    for (int l=0; l<numRelaxationMechanisms; l++) {
        update -= DT_temp * rZ[l];
    }
    update *= CbAverageZ;   
    update_temp = CaAverageZ * eZ;
    eZ = update_temp + update; 

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(eZ, hX, hY);
    }
    
    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Dviscotmem<float>;
template class KITGPI::ForwardSolver::FD2Dviscotmem<double>;
