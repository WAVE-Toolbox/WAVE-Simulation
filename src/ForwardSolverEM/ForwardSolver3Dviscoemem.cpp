#include "ForwardSolver3Dviscoemem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    return (this->estimateBoundaryMemory(config, dist, modelCoordinates, DampingBoundary, ConvPML));
}

/*! \brief Initialization of the ForwardSolver
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
void KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init3Dviscoemem");

    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefHX().getDistributionPtr() == model.getMagneticPermeability().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution */
    auto dist = wavefield.getRefHX().getDistributionPtr();

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
void KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT)
{
    scai::lama::DenseVector<ValueType> electricConductivityEffectiveOpticalAverageX;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEffectiveOpticalAverageX;
    auto const &dielectricPermittivityAverageX = model.getDielectricPermittivityAverageX();
    auto const &electricConductivityAverageX = model.getElectricConductivityAverageX();
    auto const &tauDielectricPermittivityAverageX = model.getTauDielectricPermittivityAverageX();
    auto const &tauElectricConductivityAverageX = model.getTauElectricConductivityAverageX();
    
    scai::lama::DenseVector<ValueType> electricConductivityEffectiveOpticalAverageY;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEffectiveOpticalAverageY;
    auto const &dielectricPermittivityAverageY = model.getDielectricPermittivityAverageY();
    auto const &electricConductivityAverageY = model.getElectricConductivityAverageY();
    auto const &tauDielectricPermittivityAverageY = model.getTauDielectricPermittivityAverageY();
    auto const &tauElectricConductivityAverageY = model.getTauElectricConductivityAverageY();
    
    scai::lama::DenseVector<ValueType> electricConductivityEffectiveOpticalAverageZ;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEffectiveOpticalAverageZ;
    auto const &dielectricPermittivityAverageZ = model.getDielectricPermittivityAverageZ();
    auto const &electricConductivityAverageZ = model.getElectricConductivityAverageZ();
    auto const &tauDielectricPermittivityAverageZ = model.getTauDielectricPermittivityAverageZ();
    auto const &tauElectricConductivityAverageZ = model.getTauElectricConductivityAverageZ();
    
    numRelaxationMechanisms = model.getNumRelaxationMechanisms();
    std::vector<ValueType> relaxationTime;          // = 1 / ( 2 * Pi * f_relax )
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * model.getRelaxationFrequency()[l])); // = 1 / ( 2 * Pi * f_relax )
    }
    
    dielectricPermittivityEffectiveOpticalAverageX = this->getDielectricPermittivityEffectiveOptical(dielectricPermittivityAverageX, electricConductivityAverageX, tauDielectricPermittivityAverageX, tauElectricConductivityAverageX, model.getDielectricPermittivityVacuum());
    electricConductivityEffectiveOpticalAverageX = this->getElectricConductivityEffectiveOptical(dielectricPermittivityAverageX, electricConductivityAverageX, tauDielectricPermittivityAverageX, numRelaxationMechanisms, relaxationTime);
    
    dielectricPermittivityEffectiveOpticalAverageY = this->getDielectricPermittivityEffectiveOptical(dielectricPermittivityAverageY, electricConductivityAverageY, tauDielectricPermittivityAverageY, tauElectricConductivityAverageY, model.getDielectricPermittivityVacuum());
    electricConductivityEffectiveOpticalAverageY = this->getElectricConductivityEffectiveOptical(dielectricPermittivityAverageY, electricConductivityAverageY, tauDielectricPermittivityAverageY, numRelaxationMechanisms, relaxationTime);

    dielectricPermittivityEffectiveOpticalAverageZ = this->getDielectricPermittivityEffectiveOptical(dielectricPermittivityAverageZ, electricConductivityAverageZ, tauDielectricPermittivityAverageZ, tauElectricConductivityAverageZ, model.getDielectricPermittivityVacuum());
    electricConductivityEffectiveOpticalAverageZ = this->getElectricConductivityEffectiveOptical(dielectricPermittivityAverageZ, electricConductivityAverageZ, tauDielectricPermittivityAverageZ, numRelaxationMechanisms, relaxationTime);

    CaAverageX = this->getAveragedCa(dielectricPermittivityEffectiveOpticalAverageX, electricConductivityEffectiveOpticalAverageX, DT);
    CaAverageY = this->getAveragedCa(dielectricPermittivityEffectiveOpticalAverageY, electricConductivityEffectiveOpticalAverageY, DT);
    CaAverageZ = this->getAveragedCa(dielectricPermittivityEffectiveOpticalAverageZ, electricConductivityEffectiveOpticalAverageZ, DT);
    
    CbAverageX = this->getAveragedCb(dielectricPermittivityEffectiveOpticalAverageX, electricConductivityEffectiveOpticalAverageX, DT);
    CbAverageY = this->getAveragedCb(dielectricPermittivityEffectiveOpticalAverageY, electricConductivityEffectiveOpticalAverageY, DT);
    CbAverageZ = this->getAveragedCb(dielectricPermittivityEffectiveOpticalAverageZ, electricConductivityEffectiveOpticalAverageZ, DT);
    
    Cc = this->getCc(numRelaxationMechanisms, relaxationTime, DT);
    
    CdAverageX = this->getAveragedCd(dielectricPermittivityAverageX, tauDielectricPermittivityAverageX, numRelaxationMechanisms, relaxationTime, DT);
    CdAverageY = this->getAveragedCd(dielectricPermittivityAverageY, tauDielectricPermittivityAverageY, numRelaxationMechanisms, relaxationTime, DT);
    CdAverageZ = this->getAveragedCd(dielectricPermittivityAverageZ, tauDielectricPermittivityAverageZ, numRelaxationMechanisms, relaxationTime, DT);
}

/*! \brief Running the 3-D visco-emem forward solver
 *
 * Start the 3-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timestep
 
 \par All discrete components of the 3D visco-EM wavefield in matrix-vector notation are: 
\begin{equation}
\begin{align*}
  &\vec{H}_{x}^{n+\frac{1}{2}} = \vec{H}_{x}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{myz}^{-1}\right)  \left(- \underline{D}_{\,y,f}^q \vec{E}_{z}^{n} + \underline{D}_{\,z,f}^q \vec{E}_{y}^{n} \right)\\
  &\vec{H}_{y}^{n+\frac{1}{2}} = \vec{H}_{y}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxz}^{-1}\right)  \left(- \underline{D}_{\,z,f}^q \vec{E}_{x}^{n} + \underline{D}_{\,x,f}^q \vec{E}_z^{n} \right)\\
  &\vec{H}_{z}^{n+\frac{1}{2}} = \vec{H}_{z}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxy}^{-1}\right)  \left(- \underline{D}_{\,x,f}^q \vec{E}_y^{n} + \underline{D}_{\,y,f}^q \vec{E}_x^{n} \right)\\  
  &\vec{E}_{x}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ax}\right) \vec{E}_{x}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{bx}\right) \left[ \frac{1}{\Delta h} \left( \underline{D}_{\,y,b}^q \vec{H}_z^{n+\frac{1}{2}} - \underline{D}_{\,z,b}^q \vec{H}_y^{n+\frac{1}{2}} \right) + \sum_{l=1}^L \vec{r}_{lx}^{n+\frac{1}{2}} \right]\\ 
  &\vec{E}_{y}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ay}\right) \vec{E}_{y}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{by}\right) \left[ \frac{1}{\Delta h} \left( \underline{D}_{\,z,b}^q  \vec{H}_x^{n+\frac{1}{2}} - \underline{D}_{\,x,b}^q \vec{H}_z^{n+\frac{1}{2}} \right) + \sum_{l=1}^L \vec{r}_{ly}^{n+\frac{1}{2}} \right]\\ 
  &\vec{E}_{z}^{n+1} =  \mathrm{diag} \left(\vec{C}_{az}\right) \vec{E}_{z}^{n} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{bz}\right) \left[ \frac{1}{\Delta h} \left(\underline{D}_{\,x,b}^q \vec{H}_y^{n+\frac{1}{2}} - \underline{D}_{\,y,b}^q \vec{H}_x^{n+\frac{1}{2}} \right) + \sum_{l=1}^L \vec{r}_{lz}^{n+\frac{1}{2}} \right]\\ 
  &\vec{r}_{lx}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{lx}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l x}\right) \vec{E}_{x}^{n} \\
  &\vec{r}_{ly}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{ly}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l y}\right) \vec{E}_{y}^{n}\\
  &\vec{r}_{lz}^{n+\frac{1}{2}} =  C_{c l} \vec{r}_{lz}^{n-\frac{1}{2}} + \Delta t \cdot \mathrm{diag} \left(\vec{C}_{d l z}\right) \vec{E}_{z}^{n}
\end{align*}
\end{equation}
 */

template <typename ValueType>
void KITGPI::ForwardSolver::FD3Dviscoemem<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{

    SCAI_REGION("timestep");

    /* Get references to required modelparameter */
    auto const &inverseMagneticPermeabilityAverageYZ = model.getInverseMagneticPermeabilityAverageYZ();
    auto const &inverseMagneticPermeabilityAverageXZ = model.getInverseMagneticPermeabilityAverageXZ();
    auto const &inverseMagneticPermeabilityAverageXY = model.getInverseMagneticPermeabilityAverageXY();
        
    /* Get references to required wavefields */
    auto &hX = wavefield.getRefHX();
    auto &hY = wavefield.getRefHY();
    auto &hZ = wavefield.getRefHZ();

    auto &eX = wavefield.getRefEX();
    auto &eY = wavefield.getRefEY();
    auto &eZ = wavefield.getRefEZ();
    
    auto &rX = wavefield.getRefRX();
    auto &rY = wavefield.getRefRY();
    auto &rZ = wavefield.getRefRZ();

    /* Get references to required derivatives matrices */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dzf = derivatives.getDzf();
    auto const &Dxb = derivatives.getDxb();
    auto const &Dzb = derivatives.getDzb();

    lama::Matrix<ValueType> const &DyfStaggeredX = derivatives.getDyfStaggeredX();
    lama::Matrix<ValueType> const &DybStaggeredX = derivatives.getDybStaggeredX();
    lama::Matrix<ValueType> const &DyfStaggeredZ = derivatives.getDyfStaggeredZ();
    lama::Matrix<ValueType> const &DybStaggeredZ = derivatives.getDybStaggeredZ();     

    SourceReceiverImpl::FDTD3Demem<ValueType> SourceReceiver(sources, receiver, wavefield);

    /* ----------------*/
    /* hx */
    /* ----------------*/
    update = DyfStaggeredZ * eZ;
    update_temp = Dzf * eY;
    if (useConvPML) {
        ConvPML.apply_ezy(update);
        ConvPML.apply_eyz(update_temp);
    }
    update -= update_temp;
    update *= inverseMagneticPermeabilityAverageYZ;
    hX -= update;

    /* ----------------*/
    /* hy */
    /* ----------------*/
    update = Dzf * eX;
    update_temp = Dxf * eZ;
    if (useConvPML) {
        ConvPML.apply_exz(update);
        ConvPML.apply_ezx(update_temp);
    }
    update -= update_temp;
    update *= inverseMagneticPermeabilityAverageXZ;
    hY -= update;

    /* ----------------*/
    /* hz */
    /* ----------------*/
    update = Dxf * eY;
    update_temp = DyfStaggeredX * eX;
    if (useConvPML) {
        ConvPML.apply_eyx(update);
        ConvPML.apply_exy(update_temp);
    }
    update -= update_temp;
    update *= inverseMagneticPermeabilityAverageXY;
    hZ -= update;

    /* ----------------*/
    /*  rX rY rZ       */
    /* ----------------*/
    for (int l=0; l<numRelaxationMechanisms; l++) {
        update = Cc[l] * rX[l];
        update_temp = CdAverageX[l] * eX;
        rX[l] = update_temp + update;
        
        update = Cc[l] * rY[l];
        update_temp = CdAverageY[l] * eY;
        rY[l] = update_temp + update;
        
        update = Cc[l] * rZ[l];
        update_temp = CdAverageZ[l] * eZ;
        rZ[l] = update_temp + update;
    }
    
    /* ----------------*/
    /*  ex             */
    /* ----------------*/
    update = DybStaggeredX * hZ;
    update_temp = Dzb * hY;
    if (useConvPML) {
        ConvPML.apply_hzy(update);
        ConvPML.apply_hyz(update_temp);
    }
    update -= update_temp;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        update -= DT_temp * rX[l];
    }
    update *= CbAverageX;    
    update_temp = CaAverageX * eX;
    eX = update_temp + update;
    
    /* ----------------*/
    /*  ey             */
    /* ----------------*/
    update = Dzb * hX;
    update_temp = Dxb * hZ;
    if (useConvPML) {
        ConvPML.apply_hxz(update);
        ConvPML.apply_hzx(update_temp);
    }
    update -= update_temp;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        update -= DT_temp * rY[l];
    }
    update *= CbAverageY;   
    update_temp = CaAverageY * eY;
    eY = update_temp + update;
    
    /* ----------------*/
    /*  ez             */
    /* ----------------*/
    update = Dxb * hY;
    update_temp = DybStaggeredZ * hX;
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
        DampingBoundary.apply(eZ, eY, eX, hX, hY, hZ);
    }
    
    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD3Dviscoemem<float>;
template class KITGPI::ForwardSolver::FD3Dviscoemem<double>;
