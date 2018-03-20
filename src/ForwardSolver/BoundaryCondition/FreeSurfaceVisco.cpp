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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector<ValueType> &onePlusLtauP, scai::lama::Vector<ValueType> &onePlusLtauS)
{

    lama::Vector<ValueType> const &pWaveModulus = model.getPWaveModulus();
    lama::Vector<ValueType> const &sWaveModulus = model.getSWaveModulus();
    lama::Vector<ValueType> const &tauS = model.getTauS();
    lama::Vector<ValueType> const &tauP = model.getTauP();

    /* --------------------------------------- */
    /* Apply scaling for update of Sxx and Szz */
    /* --------------------------------------- */

    auto temp = lama::eval<lama::DenseVector<ValueType>>( 2 * sWaveModulus * onePlusLtauS );
    auto temp2 = lama::eval<lama::DenseVector<ValueType>>( - pWaveModulus * onePlusLtauP );

    temp += temp2; // = ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) );

    scaleStressHorizontalUpdate = pWaveModulus;
    scaleStressHorizontalUpdate *= onePlusLtauP;           // = ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate = 1 / scaleStressHorizontalUpdate;                  // = 1 / ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate *= temp;                   // = ( ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) ) ) / ( ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate *= selectHorizontalUpdate; // set to zero everywhere besides the surface

    /* --------------------------------------- */
    /* Apply scaling for update of Rxx and Rzz */
    /* --------------------------------------- */

    temp = 2 * sWaveModulus * tauS;
    temp2 = - pWaveModulus * tauP;

    temp += temp2; // = ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP );

    scaleRelaxationHorizontalUpdate = pWaveModulus;
    scaleRelaxationHorizontalUpdate *= tauP;                   // = ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate = 1 / scaleRelaxationHorizontalUpdate;                  // = 1 / ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate *= temp;                   // = ( ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP ) ) / ( ( P-wave Modul tauP) )
    scaleRelaxationHorizontalUpdate *= selectHorizontalUpdate; // set to zero everywhere besides the surface
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfPressure(NX, NY, NZ, dist);
    derivatives.calcDyfVelocity(NX, NY, NZ, dist);
    derivatives.calcDybPressure(NX, NY, NZ, dist);
    derivatives.calcDybVelocity(NX, NY, NZ, dist);
    derivatives.DybPressure *= DT / DH;
    derivatives.DybVelocity *= DT / DH;
    derivatives.DyfPressure *= DT / DH;
    derivatives.DyfVelocity *= DT / DH;
    derivatives.Dyb.purge();
    derivatives.Dyf.purge();

    selectHorizontalUpdate.setSameValue(dist, 0.0);
    setSurfaceZero.setSameValue(dist, 1.0);

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    /* Get write access to local part of scaleHorizontalUpdate */
    auto write_selectHorizontalUpdate = hmemo::hostWriteAccess(selectHorizontalUpdate.getLocalValues());

    /* Get write access to local part of setSurfaceZero */
    auto write_setSurfaceZero = hmemo::hostWriteAccess(setSurfaceZero.getLocalValues());

    KITGPI::Acquisition::Coordinates coordinateTransformation;

    IndexType rowGlobal;
    IndexType rowLocal;

    for (IndexType i = 0; i < numLocalIndices; i++) {

        rowGlobal = read_localIndices[i];
        rowLocal = dist->global2local(rowGlobal);

        /* Determine if the current grid point is located on the surface */
        if (coordinateTransformation.locatedOnSurface(rowGlobal, NX, NY, NZ)) {

            /* Set horizontal update to +1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal] = 1.0;

            /* Set vector at surface to zero  */
            write_setSurfaceZero[rowLocal] = 0.0;
        }
    }

    HOST_PRINT(dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<double>;
