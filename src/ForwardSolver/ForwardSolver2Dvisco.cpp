#include "ForwardSolver2Dvisco.hpp"
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
void KITGPI::ForwardSolver::FD2Dvisco<ValueType>::initForwardSolver(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, const Modelparameter::Modelparameter<ValueType> &model, scai::hmemo::ContextPtr ctx, ValueType DT)
{
    /* Check if distributions of wavefields and models are the same */
    SCAI_ASSERT_ERROR(wavefield.getRefVX().getDistributionPtr()==model.getDensity().getDistributionPtr(),"Distributions of wavefields and models are not the same");
    
    /* Get distribibution */
    dmemo::DistributionPtr dist;
    lama::Vector &vX = wavefield.getRefVX();
    dist=vX.getDistributionPtr();
    
    /* Initialisation of Boundary Conditions */
    if (config.get<IndexType>("FreeSurface") && config.get<IndexType>("DampingBoundaryType")) {
	this->prepareBoundaryConditions(config, derivatives, dist, ctx);
    }
    
    /* Initialisation of auxiliary vectors and scalars*/
    common::unique_ptr<lama::Vector> tmp1(vX.newVector());  	// create new Vector(Pointer) with same configuration as vX (goes out of scope at functions end)
    updatePtr=std::move(tmp1); 					// assign tmp1 to updatePtr
    common::unique_ptr<lama::Vector> tmp2(vX.newVector()); 
    update_tempPtr=std::move(tmp2); 	
    common::unique_ptr<lama::Vector> tmp3(vX.newVector());
    vxxPtr=std::move(tmp3); 
    common::unique_ptr<lama::Vector> tmp4(vX.newVector());
    vyyPtr=std::move(tmp4);
    scai::common::unique_ptr<scai::lama::Vector> tmp5(vX.newVector());
    update2Ptr=std::move(tmp5);
    scai::common::unique_ptr<scai::lama::Vector> tmp6(vX.newVector());
    onePlusLtauPPtr=std::move(tmp6);
    scai::common::unique_ptr<scai::lama::Vector> tmp7(vX.newVector());
    onePlusLtauSPtr=std::move(tmp7);
    
    numRelaxationMechanisms = model.getNumRelaxationMechanisms();         
    relaxationTime = 1.0 / (2.0 * M_PI * model.getRelaxationFrequency());
    inverseRelaxationTime = 1.0 / relaxationTime;                         
    viscoCoeff1 = (1.0 - DT / (2.0 * relaxationTime));                    
    viscoCoeff2 = 1.0 / (1.0 + DT / (2.0 * relaxationTime));              
    DThalf = DT / 2.0;                                                    
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
void KITGPI::ForwardSolver::FD2Dvisco<ValueType>::prepareBoundaryConditions(Configuration::Configuration const &config, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx)
{

    /* Prepare Free Surface */
    if (config.get<IndexType>("FreeSurface")) {
        useFreeSurface = true;
        FreeSurface.init(dist, derivatives, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DT"), config.get<ValueType>("DH"));
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

/*! \brief Running the 2-D visco-elastic foward solver
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
 */
template <typename ValueType>
void KITGPI::ForwardSolver::FD2Dvisco<ValueType>::run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, Derivatives::Derivatives<ValueType> const &derivatives, IndexType tStart, IndexType tEnd)
{

    SCAI_REGION("timestep")

    SCAI_ASSERT_ERROR((tEnd - tStart) >= 1, " Number of time steps has to be greater than zero. ");

    /* Get references to required modelparameter */
    lama::Vector const &inverseDensity = model.getInverseDensity();
    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &sWaveModulus = model.getSWaveModulus();
    lama::Vector const &inverseDensityAverageX = model.getInverseDensityAverageX();
    lama::Vector const &inverseDensityAverageY = model.getInverseDensityAverageY();
    lama::Vector const &sWaveModulusAverageXY = model.getSWaveModulusAverageXY();
    lama::Vector const &tauSAverageXY = model.getTauSAverageXY();

    /* Get references to required wavefields */
    lama::Vector &vX = wavefield.getRefVX();
    lama::Vector &vY = wavefield.getRefVY();

    lama::Vector &Sxx = wavefield.getRefSxx();
    lama::Vector &Syy = wavefield.getRefSyy();
    lama::Vector &Sxy = wavefield.getRefSxy();

    lama::Vector &Rxx = wavefield.getRefRxx();
    lama::Vector &Ryy = wavefield.getRefRyy();
    lama::Vector &Rxy = wavefield.getRefRxy();

    /* Get references to required derivatives matrixes */
    lama::Matrix const &Dxf = derivatives.getDxf();
    lama::Matrix const &Dxb = derivatives.getDxb();

    lama::Matrix const &DybPressure = derivatives.getDybPressure();
    lama::Matrix const &DybVelocity = derivatives.getDybVelocity();
    lama::Matrix const &DyfPressure = derivatives.getDyfPressure();
    lama::Matrix const &DyfVelocity = derivatives.getDyfVelocity();

    SourceReceiverImpl::FDTD2Delastic<ValueType> SourceReceiver(sources, receiver, wavefield);

    /* Get references to auxiliary vectors */
    lama::Vector &update = *updatePtr;                          // get Reference of VectorPointer
    lama::Vector &update_temp = *update_tempPtr;                // get Reference of VectorPointer
    lama::Vector &vxx = *vxxPtr;
    lama::Vector &vyy = *vyyPtr;
    lama::Vector &update2= *update2Ptr;
    lama::Vector &onePlusLtauP=*onePlusLtauPPtr; 
    lama::Vector &onePlusLtauS=*onePlusLtauSPtr; 

    /* Get reference to required model vectors */ 
    lama::Vector const &tauS = model.getTauS();
    lama::Vector const &tauP = model.getTauP();

    onePlusLtauP = 1.0;
    onePlusLtauP += numRelaxationMechanisms * tauP;

    onePlusLtauS = 1.0;
    onePlusLtauS += numRelaxationMechanisms * tauS;

    if (useFreeSurface) {
        FreeSurface.setModelparameter(model, onePlusLtauP, onePlusLtauS);
    }

    dmemo::CommunicatorPtr comm = inverseDensity.getDistributionPtr()->getCommunicatorPtr();

    /* --------------------------------------- */
    /* Start runtime critical part             */
    /* --------------------------------------- */

    for (IndexType t = tStart; t < tEnd; t++) {

        if (t % 100 == 0 && t != 0) {
            HOST_PRINT(comm, "Calculating time step " << t << "\n");
        }

        /* ----------------*/
        /* update velocity */
        /* ----------------*/
        update = Dxf * Sxx;
        if (useConvPML) {
            ConvPML.apply_sxx_x(update);
        }

        update_temp = DybVelocity * Sxy;
        if (useConvPML) {
            ConvPML.apply_sxy_y(update_temp);
        }
        update += update_temp;
        update *= inverseDensityAverageX;
        vX += update;

        update = Dxb * Sxy;
        if (useConvPML) {
            ConvPML.apply_sxy_x(update);
        }

        update_temp = DyfVelocity * Syy;
        if (useConvPML) {
            ConvPML.apply_syy_y(update_temp);
        }
        update += update_temp;

        update *= inverseDensityAverageY;
        vY += update;

        /* ----------------*/
        /* pressure update */
        /* ----------------*/

        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        if (useConvPML) {
            ConvPML.apply_vxx(vxx);
            ConvPML.apply_vyy(vyy);
        }

        update = vxx;
        update += vyy;
        update *= pWaveModulus;

        update2 = inverseRelaxationTime * update;
        update2 *= tauP;

        Sxx += DThalf * Rxx;
        Rxx *= viscoCoeff1;
        Rxx -= update2;

        Syy += DThalf * Ryy;
        Ryy *= viscoCoeff1;
        Ryy -= update2;

        update *= onePlusLtauP;
        Sxx += update;
        Syy += update;

        /* Update Sxx and Rxx */
        vyy *= sWaveModulus;
        vyy *= 2.0;

        update2 = inverseRelaxationTime * vyy;

        update2 *= tauS;
        Rxx += update2;
        vyy *= onePlusLtauS;
        Sxx -= vyy;

        Rxx *= viscoCoeff2;
        Sxx += DThalf * Rxx;

        /* Update Syy and Ryy */
        vxx *= sWaveModulus;
        vxx *= 2.0;

        update2 = inverseRelaxationTime * vxx;
        update2 *= tauS;
        Ryy += update2;
        vxx *= onePlusLtauS;
        Syy -= vxx;

        Ryy *= viscoCoeff2;
        Syy += DThalf * Ryy;

        /* Update Sxy and Rxy*/
        Sxy += DThalf * Rxy;
        Rxy *= viscoCoeff1;

        update = DyfPressure * vX;
        if (useConvPML) {
            ConvPML.apply_vxy(update);
        }

        update_temp = Dxf * vY;
        if (useConvPML) {
            ConvPML.apply_vyx(update_temp);
        }
        update += update_temp;

        update *= sWaveModulusAverageXY;

        update2 = inverseRelaxationTime * update;
        update2 *= tauSAverageXY;
        Rxy -= update2;
        update *= onePlusLtauS;
        Sxy += update;

        Rxy *= viscoCoeff2;
        Sxy += DThalf * Rxy;

        /* Apply free surface to stress update */
        if (useFreeSurface) {
            FreeSurface.apply(vxx, update2, Sxx, Syy, Rxx, Ryy);
        }

        /* Apply the damping boundary */
        if (useDampingBoundary) {
            DampingBoundary.apply(Sxx, Syy, Sxy, vX, vY);
        }

        /* Apply source and save seismogram */
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);
    }

    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}

template class KITGPI::ForwardSolver::FD2Dvisco<float>;
template class KITGPI::ForwardSolver::FD2Dvisco<double>;
