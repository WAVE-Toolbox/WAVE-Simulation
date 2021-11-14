#include "FreeSurfaceViscoelastic.hpp"
using namespace scai;

/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 \param onePlusLtauP Parameter with ( 1 + L * tauP )
 \param onePlusLtauS Parameter with ( 1 + L * tauS )
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceViscoelastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector<ValueType> &onePlusLtauP, scai::lama::Vector<ValueType> &onePlusLtauS, ValueType DT)
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

    scaleStressVerticalUpdate = selectFreeSurface;
    scaleStressVerticalUpdate *= temp;

    scaleStressHorizontalUpdate = -1 * selectFreeSurface;
    scaleStressHorizontalUpdate *= temp;
    scaleStressHorizontalUpdate *= temp;
    scaleStressHorizontalUpdate *= temp2;

    /* --------------------------------------- */
    /* Apply scaling for update of Rxx and Rzz */
    /* --------------------------------------- */
    numRelaxationMechanisms = model.getNumRelaxationMechanisms();
    std::vector<ValueType> viscoCoeff2;
    std::vector<ValueType> relaxationTime;
    scaleRelaxationVerticalUpdate.clear();
    scaleRelaxationHorizontalUpdate.clear();
    for (int l=0; l<numRelaxationMechanisms; l++) {
        relaxationTime.push_back(1.0 / (2.0 * M_PI * model.getRelaxationFrequency()[l])); // = 1 / ( 2 * Pi * f_relax )
        viscoCoeff2.push_back(1.0 / (1.0 + DT / (2.0 * relaxationTime[l])));              // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
        scaleRelaxationVerticalUpdate.push_back(selectFreeSurface);
        scaleRelaxationHorizontalUpdate.push_back(selectFreeSurface); // set to zero everywhere besides the surface

        /* On the free surface the verical velocity derivarive can be expressed by 
        * vyy = ((2mu / pi ) -1)*(vxx+vzz)  where mu = sWaveModulus and pi = pWaveModulus 
        * tau = relaxationTime
        * The original update,
        * Rxx = (1/(1+1/(2*tau)))*(Rxx*(1-1/(2*tau))- pi*(tauP/tau)*(vxx+vyy+vzz)+mu*(tauS/tau)*(vyy+vzz))
        * will be exchanged with 
        * Rxx_new = (1/(1+1/(2*tau)))*(Rxx*(1-1/(2*tau))- pi*(tauP/tau)*(vxx+((2mu / pi ) -1)*(vxx+vzz) + vzz)+mu*(tauS/tau)*(((2mu / pi ) -1)*(vxx+vzz) + vzz))
        * The final update will be (last update has to be undone):
        * Rxx += Rxx_new - Rxx = (1/(tau+0.5))*((2mu*tauS-pi*tauP)*((2mu(1+L*tauS)/pi(1+L*tauP))-1)*(vxx+vzz) - (2mu*tauS-pi*tauP)*vyy)
        *                      = scaleRelaxationHorizontalUpdate[l]*(vxx+vzz)  - scaleRelaxationVerticalUpdate[l]*Vyy 
        * the Rzz calculation is done analogous */

        temp = 2 * sWaveModulus;
        temp2 = pWaveModulus;
        auto temp3 = lama::eval<lama::DenseVector<ValueType>>(temp * tauS);
        temp3 -= lama::eval<lama::DenseVector<ValueType>>(temp2 * tauP);

        scaleRelaxationVerticalUpdate[l] *= temp3;
        scaleRelaxationVerticalUpdate[l] *= viscoCoeff2[l];
        scaleRelaxationVerticalUpdate[l] /= relaxationTime[l];

        temp *= onePlusLtauS;
        temp2 *= onePlusLtauP;
        temp /= temp2;
        temp -= 1;

        scaleRelaxationHorizontalUpdate[l] *= temp3;
        scaleRelaxationHorizontalUpdate[l] *= temp;
        scaleRelaxationHorizontalUpdate[l] *= viscoCoeff2[l];
        scaleRelaxationHorizontalUpdate[l] /= relaxationTime[l];
    }
}

/*! \brief Initialization of the free surface
 *
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceViscoelastic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT)
{
    HOST_PRINT(dist->getCommunicatorPtr(), "", "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfFreeSurface(modelCoordinates, dist);
    derivatives.calcDybFreeSurface(modelCoordinates, dist);
    derivatives.getDybFreeSurface() *= DT;
    derivatives.getDyfFreeSurface() *= DT;

    hmemo::HArray<IndexType> ownedIndeces;
    dist->getOwnedIndexes(ownedIndeces);

    lama::VectorAssembly<ValueType> assemblyZeros;
    lama::VectorAssembly<ValueType> assemblyOnes;

    setZeroFreeSurface.setSameValue(dist, 1.0);
    selectFreeSurface.setSameValue(dist, 0.0);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {

        if (modelCoordinates.locatedOnSurface(ownedIndex)) {
            assemblyZeros.push(ownedIndex, 0);
            assemblyOnes.push(ownedIndex, 1);
        }
    }
    selectFreeSurface.fillFromAssembly(assemblyOnes);
    setZeroFreeSurface.fillFromAssembly(assemblyZeros);

    HOST_PRINT(dist->getCommunicatorPtr(), "", "Finished initializing of the free surface\n\n");
}


template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceViscoelastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceViscoelastic<double>;
