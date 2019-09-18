#include "FreeSurfaceVisco.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::~FreeSurfaceVisco(){};

/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 \param onePlusLtauP Parameter with ( 1 + L * tauP )
 \param onePlusLtauS Parameter with ( 1 + L * tauS )
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector<ValueType> &onePlusLtauP, scai::lama::Vector<ValueType> &onePlusLtauS, ValueType DT)
{
    /*This function sets scaling factors for (vxx+vzz) and vyy for the calculation of the horizontal updates without vertical derivatives */

    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();
    lama::Vector<ValueType> const &tauS = model.getTauS();
    lama::Vector<ValueType> const &tauP = model.getTauP();

    SCAI_ASSERT(sWaveModulus.min() > 0, "S wave modulus can't be zero when using image method")

    /* On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) * (vxx+vzz) where mu = sWaveModulus * (1+L*tauS) and pi = pWaveModulus * (1+L*tauP)
    * The original update,
    * sxx = pi * ( vxx+vyy+vzz ) - 2mu *(vzz+vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx + vzz - 2mu / pi * (vxx+vzz))
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx+vzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx+vzz)  - scaleVerticalUpdate*Vyy 
    * the szz calculation is done analogous */

    
    /* --------------------------------------- */
    /* Apply scaling for update of Sxx and Szz */
    /* --------------------------------------- */

    auto temp = lama::eval<lama::DenseVector<ValueType>>(-2 * sWaveModulus * onePlusLtauS);
    auto temp2 = lama::eval<lama::DenseVector<ValueType>>((pWaveModulus * onePlusLtauP));
    temp += temp2;
    temp2 = 1 / temp2;

    scaleStressVerticalUpdate = selectHorizontalUpdate;
    scaleStressVerticalUpdate *= temp;

    scaleStressHorizontalUpdate = -1 * selectHorizontalUpdate;
    scaleStressHorizontalUpdate *= temp;
    scaleStressHorizontalUpdate *= temp;
    scaleStressHorizontalUpdate *= temp2;

    /* --------------------------------------- */
    /* Apply scaling for update of Rxx and Rzz */
    /* --------------------------------------- */

    ValueType relaxationTime = 1.0 / (2.0 * M_PI * model.getRelaxationFrequency()); // = 1 / ( 2 * Pi * f_relax )
    ValueType viscoCoeff2 = 1.0 / (1.0 + DT / (2.0 * relaxationTime));              // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1 todo DT

    /* On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1)*(vxx+vzz)  where mu = sWaveModulus and pi = pWaveModulus 
    * tau = relaxationTime
    * The original update,
    * Rxx = (1/(1+1/(2*tau)))*(Rxx*(1-1/(2*tau))- pi*(tauP/tau)*(vxx+vyy+vzz)+mu*(tauS/tau)*(vyy+vzz))
    * will be exchanged with 
    * Rxx_new = (1/(1+1/(2*tau)))*(Rxx*(1-1/(2*tau))- pi*(tauP/tau)*(vxx+((2mu / pi ) -1)*(vxx+vzz) + vzz)+mu*(tauS/tau)*(((2mu / pi ) -1)*(vxx+vzz) + vzz))
    * The final update will be (last update has to be undone):
    * Rxx += Rxx_new - Rxx = (1/(tau+0.5))*((2mu*tauS-pi*tauP)*((2mu(1+L*tauS)/pi(1+L*tauP))-1)*(vxx+vzz) - (2mu*tauS-pi*tauP)*vyy)
    *                      = scaleRelaxationHorizontalUpdate*(vxx+vzz)  - scaleRelaxationVerticalUpdate*Vyy 
    * the Rzz calculation is done analogous */

    temp = 2 * sWaveModulus;
    temp2 = pWaveModulus;
    auto temp3 = lama::eval<lama::DenseVector<ValueType>>(temp * tauS);
    temp3 -= lama::eval<lama::DenseVector<ValueType>>(temp2 * tauP);

    scaleRelaxationVerticalUpdate = selectHorizontalUpdate;
    scaleRelaxationVerticalUpdate *= temp3;
    scaleRelaxationVerticalUpdate *= viscoCoeff2;
    scaleRelaxationVerticalUpdate /= relaxationTime;

    temp *= onePlusLtauS;
    temp2 *= onePlusLtauP;
    temp /= temp2;
    temp -= 1;

    scaleRelaxationHorizontalUpdate = selectHorizontalUpdate; // set to zero everywhere besides the surface
    scaleRelaxationHorizontalUpdate *= temp3;
    scaleRelaxationHorizontalUpdate *= temp;
    scaleRelaxationHorizontalUpdate *= viscoCoeff2;
    scaleRelaxationHorizontalUpdate /= relaxationTime;
}

/*! \brief Initialitation of the free surface
 *
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "", "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfFreeSurface(modelCoordinates, dist);
    derivatives.calcDybFreeSurface(modelCoordinates, dist);
    derivatives.getDybFreeSurface() *= DT;
    derivatives.getDyfFreeSurface() *= DT;

    selectHorizontalUpdate.setSameValue(dist, 0.0);
    setSurfaceZero.setSameValue(dist, 1.0);

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    scai::lama::DenseVector<ValueType> temp2(dist, 0.0);
    /* Get write access to local part of temp2 */
    auto write_selectHorizontalUpdate = hmemo::hostWriteAccess(temp2.getLocalValues());

    scai::lama::DenseVector<ValueType> temp(dist, 1.0);

    /* Get write access to local part of temp */
    auto write_setSurfaceZero = hmemo::hostWriteAccess(temp.getLocalValues());

    IndexType rowGlobal;
    IndexType rowLocal;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        rowGlobal = read_localIndices[i];
        rowLocal = dist->global2Local(rowGlobal);

        /* Determine if the current grid point is located on the surface */
        if (modelCoordinates.locatedOnSurface(rowGlobal)) {

            /* Set horizontal update to +1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal] = 1.0;

            /* Set vector at surface to zero  */
            write_setSurfaceZero[rowLocal] = 0.0;
        }
    }
    write_setSurfaceZero.release();
    write_selectHorizontalUpdate.release();
    setSurfaceZero = temp;
    selectHorizontalUpdate = temp2;

    HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished initializing of the free surface\n\n");
}

/*! \brief Set Ryy at the free surface to zero
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param Ryy Ryy wavefield (relaxation)
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::setMemoryVariableToZero(scai::lama::Vector<ValueType> &Ryy)
{

    SCAI_ASSERT_DEBUG(active, " FreeSurface is not active ");

    Ryy *= setSurfaceZero; // Set at the free surface to zero
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<double>;
