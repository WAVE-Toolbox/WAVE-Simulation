#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "../Derivatives/FDTD3D.hpp"
#include "FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D visco free surface
            template <typename ValueType>
            class FreeSurfaceVisco : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceVisco(){};

                //! Default destructor
                virtual ~FreeSurfaceVisco() = 0;

                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, lama::Vector &onePlusLtauP, lama::Vector &onePlusLtauS);

              protected:
                using FreeSurface<ValueType>::active;

                lama::DenseVector<ValueType> setSurfaceZero;         //!< Vector, which sets the wavefields at the surface to zero
                lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                lama::DenseVector<ValueType> scaleStressHorizontalUpdate;     //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the stress update
                lama::DenseVector<ValueType> scaleRelaxationHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the update of the relaxation mechanism
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */

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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, lama::Vector &onePlusLtauP, lama::Vector &onePlusLtauS)
{

    lama::Vector const &pWaveModulus = model.getPWaveModulus();
    lama::Vector const &sWaveModulus = model.getSWaveModulus();
    lama::Vector const &tauS = model.getTauS();
    lama::Vector const &tauP = model.getTauP();

    lama::DenseVector<ValueType> temp(sWaveModulus.getDistributionPtr());
    lama::DenseVector<ValueType> temp2(sWaveModulus.getDistributionPtr());

    /* --------------------------------------- */
    /* Apply scaling for update of Sxx and Szz */
    /* --------------------------------------- */

    temp = 2 * sWaveModulus;
    temp.scale(onePlusLtauS);

    temp2 = -1.0 * pWaveModulus;
    temp2.scale(onePlusLtauP);

    temp += temp2; // = ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) );

    scaleStressHorizontalUpdate = pWaveModulus;
    scaleStressHorizontalUpdate.scale(onePlusLtauP);           // = ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.invert();                      // = 1 / ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.scale(temp);                   // = ( ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) ) ) / ( ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.scale(selectHorizontalUpdate); // set to zero everywhere besides the surface

    /* --------------------------------------- */
    /* Apply scaling for update of Rxx and Rzz */
    /* --------------------------------------- */

    temp = 2 * sWaveModulus;
    temp.scale(tauS);

    temp2 = -1.0 * pWaveModulus;
    temp2.scale(tauP);

    temp += temp2; // = ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP );

    scaleRelaxationHorizontalUpdate = pWaveModulus;
    scaleRelaxationHorizontalUpdate.scale(tauP);                   // = ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate.invert();                      // = 1 / ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate.scale(temp);                   // = ( ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP ) ) / ( ( P-wave Modul tauP) )
    scaleRelaxationHorizontalUpdate.scale(selectHorizontalUpdate); // set to zero everywhere besides the surface
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH)
{

    HOST_PRINT(dist->getCommunicatorPtr(), "Initialization of the free surface...\n");

    active = true;

    derivatives.useFreeSurface = true;
    derivatives.calcDyfPressure(NX, NY, NZ, dist);
    derivatives.calcDyfVelocity(NX, NY, NZ, dist);
    derivatives.calcDybPressure(NX, NY, NZ, dist);
    derivatives.calcDybVelocity(NX, NY, NZ, dist);
    derivatives.DybPressure.scale(lama::Scalar(DT / DH));
    derivatives.DybVelocity.scale(lama::Scalar(DT / DH));
    derivatives.DyfPressure.scale(lama::Scalar(DT / DH));
    derivatives.DyfVelocity.scale(lama::Scalar(DT / DH));
    derivatives.Dyb.purge();
    derivatives.Dyf.purge();

    selectHorizontalUpdate.allocate(dist);
    selectHorizontalUpdate = 0.0;

    setSurfaceZero.allocate(dist);
    setSurfaceZero = 1.0;

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);                          /* get local indices based on used distribution */
    IndexType numLocalIndices = localIndices.size();              // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices

    /* Get write access to local part of scaleHorizontalUpdate */
    utilskernel::LArray<ValueType> *selectHorizontalUpdate_LA = &selectHorizontalUpdate.getLocalValues();
    hmemo::WriteAccess<ValueType> write_selectHorizontalUpdate(*selectHorizontalUpdate_LA);

    /* Get write access to local part of setSurfaceZero */
    utilskernel::LArray<ValueType> *setSurfaceZero_LA = &setSurfaceZero.getLocalValues();
    hmemo::WriteAccess<ValueType> write_setSurfaceZero(*setSurfaceZero_LA);

    KITGPI::Acquisition::Coordinates<ValueType> coordinateTransformation;

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
    read_localIndices.release();
    write_selectHorizontalUpdate.release();
    write_setSurfaceZero.release();
    HOST_PRINT(dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n");
}
