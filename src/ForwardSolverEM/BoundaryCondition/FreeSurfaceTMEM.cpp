#include "FreeSurfaceTMEM.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceTMEM<ValueType>::~FreeSurfaceTMEM(){};

/*! \brief Scale horizontal  and vertical updates with model parameter
 * this will be used to exchange vertical with horizontal derivatives for the horizontal updates on the free surface 
 *
 *
 \param model which is used during forward modelling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceTMEM<ValueType>::setModelparameter(Modelparameter::ModelparameterEM<ValueType> const &model)
{
    /*This function sets scaling factors for (vxx+hzz) and vyy for the calculation of the horizontal updates without vertical derivatives 
    * On the free surface the verical velocity derivarive can be expressed by 
    * vyy = ((2mu / pi ) -1) (vxx+hzz) where mu = dielectricPermittivityEM and pi = velocivityEM
    * The original update,
    * sxx = pi * ( vxx+vyy+hzz ) - 2mu *( hzz+vyy )
    * will be exchanged with 
    * sxx_new = 2mu * (2vxx + hzz - 2mu / pi * (vxx+hzz))
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx + hzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx + hzz)  -  scaleVerticalUpdate*Hyy 
    * The update of szz is calculated the same way   
    */
/*
    lama::Vector<ValueType> const &velocivityEM = model.getVelocityEM();
    lama::Vector<ValueType> const &dielectricPermittivityEM = model.getDielectricPermittivityEM();

    SCAI_ASSERT(dielectricPermittivityEM.min() > 0, "S wave modulus can't be zero when using image method")

    auto temp = lama::eval<lama::DenseVector<ValueType>>(velocivityEM - 2 * dielectricPermittivityEM);

    scaleVerticalUpdate = selectFreeSurface;
    scaleVerticalUpdate *= temp;

    temp *= temp;
    temp /= velocivityEM;
    temp *= -1;

    scaleHorizontalUpdate = selectFreeSurface;

    scaleHorizontalUpdate *= temp;*/
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceTMEM<ValueType>::init(scai::dmemo::DistributionPtr dist, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType /*DT*/)
{
//     dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();
// 
//     HOST_PRINT(comm, "", "Initialization of the free surface...\n");
// 
//     active = true;
// 
//     derivatives.useFreeSurface = true;
// 
//     hmemo::HArray<IndexType> ownedIndeces;
//     dist->getOwnedIndexes(ownedIndeces);
// 
//     lama::VectorAssembly<ValueType> assemblyZeros;
//     lama::VectorAssembly<ValueType> assemblyOnes;
// 
//     setZeroFreeSurface.setSameValue(dist, 1.0);
//     selectFreeSurface.setSameValue(dist, 0.0);
// 
//     for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndeces)) {
// 
//         if (modelCoordinates.locatedOnSurface(ownedIndex)) {
//             assemblyZeros.push(ownedIndex, 0);
//             assemblyOnes.push(ownedIndex, 1);
//         }
//     }
//     selectFreeSurface.fillFromAssembly(assemblyOnes);
//     setZeroFreeSurface.fillFromAssembly(assemblyZeros);
// 
//     HOST_PRINT(comm, "", "Finished initializing of the free surface\n\n");
}


template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceTMEM<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceTMEM<double>;
