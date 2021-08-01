#include "ForwardSolver3Demem.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD3Demem<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Demem<ValueType>::initForwardSolver(Configuration::Configuration const &config, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Wavefields::WavefieldsEM<ValueType> &wavefield, Modelparameter::ModelparameterEM<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
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
void KITGPI::ForwardSolver::FD3Demem<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Demem<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD3Demem<ValueType>::prepareForModelling(Modelparameter::ModelparameterEM<ValueType> const &model, ValueType DT)
{
    auto const &electricConductivityAverageX = model.getElectricConductivityAverageX();
    auto const &electricConductivityAverageY = model.getElectricConductivityAverageY();
    auto const &electricConductivityAverageZ = model.getElectricConductivityAverageZ();
    auto const &dielectricPermittivityAverageX = model.getDielectricPermittivityAverageX();
    auto const &dielectricPermittivityAverageY = model.getDielectricPermittivityAverageY();
    auto const &dielectricPermittivityAverageZ = model.getDielectricPermittivityAverageZ();

    CaAverageX = this->getAveragedCa(dielectricPermittivityAverageX, electricConductivityAverageX, DT);
    CaAverageY = this->getAveragedCa(dielectricPermittivityAverageY, electricConductivityAverageY, DT);
    CaAverageZ = this->getAveragedCa(dielectricPermittivityAverageZ, electricConductivityAverageZ, DT);
    
    CbAverageX = this->getAveragedCb(dielectricPermittivityAverageX, electricConductivityAverageX, DT);
    CbAverageY = this->getAveragedCb(dielectricPermittivityAverageY, electricConductivityAverageY, DT);
    CbAverageZ = this->getAveragedCb(dielectricPermittivityAverageZ, electricConductivityAverageZ, DT);
}

/*! \brief Running the 3-D emem foward solver
 *
 * Start the 3-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timestep
 *
 \par All discrete components of the 3D EM wavefield in matrix-vector notation are: 
\begin{equation}
\begin{align*}
  &\vec{H}_{x}^{n+\frac{1}{2}} = \vec{H}_{x}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{myz}^{-1}\right)  \left(- \underline{D}_{\,y,f}^q \vec{E}_{z}^{n} + \underline{D}_{\,z,f}^q \vec{E}_{y}^{n} \right)\\
  &\vec{H}_{y}^{n+\frac{1}{2}} = \vec{H}_{y}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxz}^{-1}\right)  \left(- \underline{D}_{\,z,f}^q \vec{E}_{x}^{n} + \underline{D}_{\,x,f}^q \vec{E}_z^{n} \right)\\
  &\vec{H}_{z}^{n+\frac{1}{2}} = \vec{H}_{z}^{n-\frac{1}{2}} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{\mu}_{mxy}^{-1}\right)  \left(- \underline{D}_{\,x,f}^q \vec{E}_y^{n} + \underline{D}_{\,y,f}^q \vec{E}_x^{n} \right)\\
  &\vec{E}_{x}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ax}\right) \vec{E}_{x}^{n} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{bx}\right) \left( \underline{D}_{\,y,b}^q \vec{H}_z^{n+\frac{1}{2}} - \underline{D}_{\,z,b}^q \vec{H}_y^{n+\frac{1}{2}} \right)\\
  &\vec{E}_{y}^{n+1} =  \mathrm{diag} \left(\vec{C}_{ay}\right) \vec{E}_{y}^{n} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{by}\right) \left( \underline{D}_{\,z,b}^q  \vec{H}_x^{n+\frac{1}{2}} - \underline{D}_{\,x,b}^q \vec{H}_z^{n+\frac{1}{2}} \right)\\
  &\vec{E}_{z}^{n+1} =  \mathrm{diag} \left(\vec{C}_{az}\right) \vec{E}_{z}^{n} + \frac{\Delta t}{\Delta h} \cdot  \mathrm{diag} \left(\vec{C}_{bz}\right) \left(\underline{D}_{\,x,b}^q \vec{H}_y^{n+\frac{1}{2}} - \underline{D}_{\,y,b}^q \vec{H}_x^{n+\frac{1}{2}} \right)
\end{align*}
\end{equation}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD3Demem<ValueType>::run(Acquisition::AcquisitionGeometryEM<ValueType> &receiver, Acquisition::AcquisitionGeometryEM<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
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

    /* Get references to required derivatives matrixes */
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
    /*  ex */
    /* ----------------*/
    update = DybStaggeredX * hZ;
    update_temp = Dzb * hY;
    if (useConvPML) {
        ConvPML.apply_hzy(update);
        ConvPML.apply_hyz(update_temp);
    }
    update -= update_temp;
    update *= CbAverageX;    
    update_temp = CaAverageX * eX;
    eX = update_temp + update;
    
    /* ----------------*/
    /*  ey */
    /* ----------------*/
    update = Dzb * hX;
    update_temp = Dxb * hZ;
    if (useConvPML) {
        ConvPML.apply_hxz(update);
        ConvPML.apply_hzx(update_temp);
    }
    update -= update_temp;
    update *= CbAverageY;   
    update_temp = CaAverageY * eY;
    eY = update_temp + update;
    
    /* ----------------*/
    /*  ez */
    /* ----------------*/
    update = Dxb * hY;
    update_temp = DybStaggeredZ * hX;
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
        DampingBoundary.apply(eZ, eY, eX, hX, hY, hZ);
    }

    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template <typename ValueType>
void KITGPI::ForwardSolver::FD3Demem<ValueType>::runAdjoint(Acquisition::AcquisitionGeometryEM<ValueType> &receiver, Acquisition::AcquisitionGeometryEM<ValueType> const &sources, Modelparameter::ModelparameterEM<ValueType> const &model, Wavefields::WavefieldsEM<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{

}

template class KITGPI::ForwardSolver::FD3Demem<float>;
template class KITGPI::ForwardSolver::FD3Demem<double>;
