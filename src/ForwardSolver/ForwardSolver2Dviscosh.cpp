#include "ForwardSolver2Dviscosh.hpp"
using namespace scai;

template <typename ValueType>
ValueType KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
void KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    SCAI_REGION("ForwardSolver.init2Dviscosh");
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVZ().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    auto dist = wavefield.getRefVZ().getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, modelCoordinates, derivatives, dist, ctx);
    }

    /* allocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);
    update2.allocate(dist);
    onePlusLtauS.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
    update2.setContextPtr(ctx);
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
void KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{
    this->prepareBoundaries(config, modelCoordinates, derivatives, dist, ctx, FreeSurface, DampingBoundary, ConvPML);
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType /*DT*/)
{
    /* Get reference to required model vectors */
    lama::Vector<ValueType> const &tauS = model.getTauS();

    onePlusLtauS = 1.0;
    onePlusLtauS += numRelaxationMechanisms * tauS;
}

/*! \brief Running the 2-D viscosh foward solver
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
 * The update equations for velocity, \f$v_i\f$, and stress, \f$\sigma_{ij}\f$, are implemented as follows where \f$\mu\f$ is the S-wave modulus and \f$\rho_{inv}\f$ the inverse density. Note that the scaling with the temporal and spatial discretization is included in the derivative matrices. The velocity update is executed first followed by the stress update and finally the source term is added. If a free surface is chosen, the derivative matrices will be adapted to satisfy the free surface condition.
 *
 \f{eqnarray*}{
	\vec{v}_z &+=& \frac{\Delta t}{\Delta h} ~ \mathrm{diag} \left( \vec{\rho}_\mathrm{inv}^{\,T} \right) \cdot \left( \underline{D}_{\,x,b}^q \vec{\sigma}_{xz} + \underline{D}_{\,y,b}^q \vec{\sigma}_{yz} \right)\\
 \f}
 \f{eqnarray*}{
	\vec{\sigma}_{xz} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,x,f}^q \vec{v}_z \right)\\
	\vec{\sigma}_{yz} &+=& \frac{\Delta t}{\Delta h}~ \mathrm{diag} \left( \vec{\mu}^{\,T} \right) \cdot \left( \underline{D}_{\,y,f}^q \vec{v}_z \right)
 \f}
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dviscosh<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t)
{
    SCAI_REGION("ForwardSolver.timestep2Dviscosh");

    /* Get references to required modelparameter */
    auto const &inverseDensity = model.getInverseDensity();
    auto const &sWaveModulusAverageXZ = model.getSWaveModulusAverageXZ();
    auto const &sWaveModulusAverageYZ = model.getSWaveModulusAverageYZ();
    auto const &tauSAverageXZ = model.getTauSAverageXZ();
    auto const &tauSAverageYZ = model.getTauSAverageYZ();

    /* Get references to required wavefields */
    auto &vZ = wavefield.getRefVZ();
    auto &Sxz = wavefield.getRefSxz();
    auto &Syz = wavefield.getRefSyz();
    auto &Ryz = wavefield.getRefRyz();
    auto &Rxz = wavefield.getRefRxz();
//     std::cout<< "Rxz[l] = " << Rxz[l] <<std::endl;
//     std::cout<< "viscoCoeff1[l] = " << viscoCoeff1[l] <<std::endl;

    /* Get references to required derivatives matrices */
    auto const &Dxf = derivatives.getDxf();
    auto const &Dxb = derivatives.getDxb();
    auto const &Dyf = derivatives.getDyf();
    auto const &Dyb = derivatives.getDyb();
    
    auto const &DybFreeSurface = derivatives.getDybFreeSurface();
        
    SourceReceiverImpl::FDTD2Dsh<ValueType> SourceReceiver(sources, receiver, wavefield);
//     if (useFreeSurface) {
//         SCAI_ASSERT(useFreeSurface != true, " Stress-image method is not implemented for Love-Waves ");
//         //         FreeSurface.setModelparameter(model);
//     }

    /* ----------------*/
    /* update velocity */
    /* ----------------*/

    update = Dxb * Sxz;

    if (useFreeSurface == 1) {
        /* Apply image method */
        update_temp = DybFreeSurface * Syz;
    } else {
        update_temp = Dyb * Syz;
    }

    if (useConvPML) {
        ConvPML.apply_sxz_x(update);
        ConvPML.apply_syz_y(update_temp);
    }
    update += update_temp;

    update *= inverseDensity;
    vZ += update;

    /* ----------------*/
    /*  update stress  */
    /* ----------------*/
    
    /* Update Sxz and Rxz */
    update = Dxf * vZ;
    if (useConvPML) {
        ConvPML.apply_vzx(update);
    }
    update *= sWaveModulusAverageXZ;
    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Sxz += DThalf * Rxz[l];
        Rxz[l] *= viscoCoeff1[l];
        
        update2 = inverseRelaxationTime[l] * update;
        update2 *= tauSAverageXZ;
        Rxz[l] -= update2;

        Rxz[l] *= viscoCoeff2[l];
        Sxz += DThalf * Rxz[l];
    }
    update *= onePlusLtauS;
    Sxz += update;

    /* Update Syz and Ryz */
    update = Dyf * vZ;
    if (useConvPML) {
        ConvPML.apply_vzy(update);
    }
    update *= sWaveModulusAverageYZ;
    
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Syz += DThalf * Ryz[l];
        Ryz[l] *= viscoCoeff1[l];
        
        update2 = inverseRelaxationTime[l] * update;
        update2 *= tauSAverageYZ;
        Ryz[l] -= update2;

        Ryz[l] *= viscoCoeff2[l];
        Syz += DThalf * Ryz[l];
    }
    update *= onePlusLtauS;
    Syz += update;
    /* Apply free surface to stress update */
//     if (useFreeSurface) {
//         SCAI_ASSERT(useFreeSurface != true, " Stress-image method is not implemented for Love-Waves ");
//         //       FreeSurface.apply(vZ, Sxz, Syz);
//     }

    /* Apply the damping boundary */
    if (useDampingBoundary) {
        DampingBoundary.apply(Sxz, Syz, vZ);
    }

    /* Apply source and save seismogram */
    SourceReceiver.applySource(t);
    SourceReceiver.gatherSeismogram(t);
}

template class KITGPI::ForwardSolver::FD2Dviscosh<float>;
template class KITGPI::ForwardSolver::FD2Dviscosh<double>;
