#include "FreeSurfaceElastic.hpp"
using namespace scai;

/*! \brief Scale horizontal  and vertical updates with model parameter
 * this will be used to exchange vertical with horizontal derivatives for the horizontal updates on the free surface 
 *
 *
 \param model which is used during forward modelling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model)
{
    /*This function sets scaling factors for (vxx+vzz) and vyy for the calculation of the horizontal updates without vertical derivatives 
    * On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) (vxx+vzz) where mu = sWaveModulus and pi = pWaveModulus
    * The original update,
    * sxx = pi * ( vxx+vyy+vzz ) - 2mu *( vzz+vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx + vzz - 2mu / pi * (vxx+vzz))
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx + vzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx + vzz)  -  scaleVerticalUpdate*Vyy 
    * The update of szz is calculated the same way   
    */

    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();

    if (sWaveModulus.min() <= 0) {
        std::cout << "sWaveModulus.min() = " << sWaveModulus.min() << std::endl;
    }
    
    SCAI_ASSERT(sWaveModulus.min() > 0, "S wave modulus can't be zero when using image method")

    auto temp = lama::eval<lama::DenseVector<ValueType>>(pWaveModulus - 2 * sWaveModulus);

    scaleVerticalUpdate = selectFreeSurface;
    scaleVerticalUpdate *= temp;

    temp *= temp;
    temp /= pWaveModulus;
    temp *= -1;

    scaleHorizontalUpdate = selectFreeSurface;

    scaleHorizontalUpdate *= temp;
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType /*DT*/)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;

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

    HOST_PRINT(comm, "", "Finished initializing of the free surface\n\n");
}


template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<double>;
