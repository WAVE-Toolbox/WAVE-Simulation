#include "FreeSurfaceElastic.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::~FreeSurfaceElastic(){};

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
    * sxx = pi * ( vxx+vyy+vzz ) - 2mu *vyy 
    * will be exchanged with 
    * sxx_new = 2mu * (2 - 2mu / pi )* (vxx+vzz)
    * The final update will be (last update has to be undone):
    * sxx += sxx_new - sxx = -(pi-2mu)*(pi-2mu)/pi*(vxx+vzz)) - (pi-2mu)*vyy
    *                      = scaleHorizontalUpdate*(vxx + vzz)  - scaleVerticalUpdate*Vyy */

    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();

    SCAI_ASSERT(sWaveModulus.min() > 0, "S wave modulus can't be zero when using image method")

    auto temp = lama::eval<lama::DenseVector<ValueType>>(pWaveModulus - 2 * sWaveModulus);

    scaleVerticalUpdate = selectHorizontalUpdate;
    scaleVerticalUpdate *= temp;

    temp *= temp;
    temp /= pWaveModulus;
    temp *= -1;

    scaleHorizontalUpdate = selectHorizontalUpdate;

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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfFreeSurface(modelCoordinates, dist);
    //     derivatives.calcDybFreeSurface(modelCoordinates, dist);
    derivatives.calcDybStaggeredXFreeSurface(modelCoordinates, dist);
    derivatives.calcDybStaggeredZFreeSurface(modelCoordinates, dist);

    derivatives.DyfFreeSurface *= DT;
    //     derivatives.DybFreeSurface *= DT;
    derivatives.DybStaggeredXFreeSurface *= DT;
    derivatives.DybStaggeredZFreeSurface *= DT;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    lama::DenseVector<ValueType> temp(dist, 0.0);
    /* Get write access to local part of scaleHorizontalUpdate */
    auto write_selectHorizontalUpdate = hostWriteAccess(temp.getLocalValues());

    IndexType rowGlobal;
    IndexType rowLocal;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        rowGlobal = read_localIndices[i];
        rowLocal = dist->global2Local(rowGlobal);

        /* Determine if the current grid point is located on the surface */
        if (modelCoordinates.locatedOnSurface(rowGlobal)) {

            /* Set horizontal update to 1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal] = 1.0;
        }
    }
    write_selectHorizontalUpdate.release();
    selectHorizontalUpdate = temp;

    HOST_PRINT(comm, "", "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<double>;
