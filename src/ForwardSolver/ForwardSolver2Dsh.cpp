#include "ForwardSolver2Dsh.hpp"
using namespace scai;

/*! \brief Initialitation of the ForwardSolver
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param wavefield Wavefields for the modelling
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dsh<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, scai::hmemo::ContextPtr ctx, ValueType /*DT*/)
{
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVZ().getDistributionPtr() == model.getDensity().getDistributionPtr(), "Distributions of wavefields and models are not the same");

    /* Get distribibution of the wavefields */
    dmemo::DistributionPtr dist;
    lama::Vector<ValueType> &vZ = wavefield.getRefVZ();
    dist = vZ.getDistributionPtr();

    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") || config.get<IndexType>("DampingBoundary")) {
        this->prepareBoundaryConditions(config, derivatives, dist, ctx);
    }

    /* aalocation of auxiliary vectors*/
    update.allocate(dist);
    update_temp.allocate(dist);

    update.setContextPtr(ctx);
    update_temp.setContextPtr(ctx);
}

/*! \brief Initialitation of the boundary conditions
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dsh<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> & /*derivatives*/, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{

      useFreeSurface = config.get<IndexType>("FreeSurface");
  
    /* Prepare Free Surface */
    if (useFreeSurface==1) {
        useFreeSurface = true;
        SCAI_ASSERT(useFreeSurface != true, " Image method is not implemented for Love-Waves ");
           //    FreeSurface.init(dist, derivatives, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<ValueType>("DH"));
    }

    /* Prepare Damping Boundary */
    if (config.get<IndexType>("DampingBoundary") == 1) {
        if (config.get<IndexType>("DampingBoundaryType") == 1) {
            useDampingBoundary = true;
            DampingBoundary.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"), useFreeSurface);
        }
        if (config.get<IndexType>("DampingBoundaryType") == 2) {
            useConvPML = true;
            ConvPML.init(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<IndexType>("DH"), config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("NPower"), config.get<ValueType>("KMaxCPML"), config.get<ValueType>("CenterFrequencyCPML"), config.get<ValueType>("VMaxCPML"), useFreeSurface);
        }
    }
}

/*! \brief resets PML (use after each modelling!)
 *
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dsh<ValueType>::resetCPML()
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
void KITGPI::ForwardSolver::FD2Dsh<ValueType>::prepareForModelling(Modelparameter::Modelparameter<ValueType> const &/*model*/, ValueType /*DT*/)
{
    if (useFreeSurface==1) {
              SCAI_ASSERT(useFreeSurface != true, " Image method is not implemented for Love-Waves ");
    //    FreeSurface.setModelparameter(model);
    }
}

/*! \brief Running the 2-D sh foward solver
 *
 * Start the 2-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param t current timestep
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dsh<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType t)
{

    SCAI_REGION("timestep")


    /* Get references to required modelparameter */
    lama::Vector<ValueType> const &inverseDensity = model.getInverseDensity();
    lama::Vector<ValueType> const &sWaveModulusAverageXZ = model.getSWaveModulusAverageXZ();
    lama::Vector<ValueType> const &sWaveModulusAverageYZ = model.getSWaveModulusAverageYZ();

    /* Get references to required wavefields */
    lama::Vector<ValueType> &vZ = wavefield.getRefVZ();
    lama::Vector<ValueType> &Sxz = wavefield.getRefSxz();
    lama::Vector<ValueType> &Syz = wavefield.getRefSyz();

    /* Get references to required derivatives matrixes */
    lama::Matrix<ValueType> const &Dxf = derivatives.getDxf();
    lama::Matrix<ValueType> const &Dxb = derivatives.getDxb();
    lama::Matrix<ValueType> const &Dyb = derivatives.getDyb();
    lama::Matrix<ValueType> const &Dyf = derivatives.getDyf();
 //   lama::Matrix const &DybVelocity = derivatives.getDybVelocity();
    
    SourceReceiverImpl::FDTD2Dsh<ValueType> SourceReceiver(sources, receiver, wavefield);


    if (useFreeSurface) {
        SCAI_ASSERT(useFreeSurface != true, " Image method is not implemented for Love-Waves ");
   //         FreeSurface.setModelparameter(model);
    }

        /* ----------------*/
        /* update velocity */
        /* ----------------*/

        update = Dxb * Sxz;
        update_temp = Dyb * Syz;
	
        if (useConvPML) {
	    ConvPML.apply_vxx(update);
            ConvPML.apply_vyy(update_temp);
        }
        
        update += update_temp;

        update *= inverseDensity;
        vZ += update;

        /* ----------------*/
        /* pressure update */
        /* ----------------*/

        update = Dxf*vZ;
        if (useConvPML) {
            ConvPML.apply_sxx_x(update);
        }
	
        update *= sWaveModulusAverageXZ;

        Sxz += update;

        update = Dyf*vZ;
        if (useConvPML) {
            ConvPML.apply_syy_y(update);
        }
        update *= sWaveModulusAverageYZ;

        Syz += update;

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            SCAI_ASSERT(useFreeSurface != true, " Image method is not implemented for Love-Waves ");
                //       FreeSurface.apply(vZ, Sxz, Syz);
        }

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(Sxz, Syz, vZ);
        }

        /* Apply source and save seismogram */
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);

}

template class KITGPI::ForwardSolver::FD2Dsh<float>;
template class KITGPI::ForwardSolver::FD2Dsh<double>;
