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
    
    SCAI_ASSERT(sWaveModulus.min()>0,"S wave modulus can't be zero when using image method")

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
 \param NX Number of grid points in X-Direction
 \param NY Number of grid points in Y-Direction (Depth)
 \param NZ Number of grid points in Z-Direction
 \param DT Temporal Sampling
 \param DH Distance between grid points
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfFreeSurface(NX, NY, NZ, dist);
    derivatives.calcDybFreeSurface(NX, NY, NZ, dist);

    derivatives.DybFreeSurface *= DT / DH;
    derivatives.DyfFreeSurface *= DT / DH;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    lama::DenseVector<ValueType> temp(dist, 0.0);
    /* Get write access to local part of scaleHorizontalUpdate */
    auto write_selectHorizontalUpdate = hostWriteAccess(temp.getLocalValues());

    KITGPI::Acquisition::Coordinates coordinateTransformation;

    IndexType rowGlobal;
    IndexType rowLocal;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        rowGlobal = read_localIndices[i];
        rowLocal = dist->global2local(rowGlobal);

        /* Determine if the current grid point is located on the surface */
        if (coordinateTransformation.locatedOnSurface(rowGlobal, NX, NY, NZ)) {

            /* Set horizontal update to 1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal] = 1.0;
        }
    }
    write_selectHorizontalUpdate.release();
    selectHorizontalUpdate = temp;

    HOST_PRINT(comm, "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<double>;
