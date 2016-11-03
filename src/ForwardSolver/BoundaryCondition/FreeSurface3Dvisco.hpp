#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../Derivatives/FD3D.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurface.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D visco free surface
            template<typename ValueType>
            class FreeSurface3Dvisco : public FreeSurface<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface3Dvisco(){};
                
                //! Default destructor
                ~FreeSurface3Dvisco(){};
                
                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH);
                
                void setModelparameter(Modelparameter::Modelparameter<ValueType>& model,lama::DenseVector<ValueType>& tauP,lama::DenseVector<ValueType>& tauS,lama::DenseVector<ValueType>& onePlusLtauP,lama::DenseVector<ValueType>& onePlusLtauS);
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::Vector& temp, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz, lama::DenseVector<ValueType>& Rxx, lama::DenseVector<ValueType>& Ryy, lama::DenseVector<ValueType>& Rzz);
                
            private:
                
                using FreeSurface<ValueType>::active;
                
                lama::DenseVector<ValueType> setSurfaceZero; //!< Vector, which sets the wavefields at the surface to zero
                lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                lama::DenseVector<ValueType> scaleStressHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the stress update
                lama::DenseVector<ValueType> scaleRelaxationHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the update of the relaxation mechanism

            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */



/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 \param tauP Tau-parameter for P-wave modulus
 \param tauS Tau-parameter for S-wave modulus
 \param onePlusLtauP Parameter with ( 1 + L * tauP )
 \param onePlusLtauS Parameter with ( 1 + L * tauS )
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType>& model,lama::DenseVector<ValueType>& tauP,lama::DenseVector<ValueType>& tauS,lama::DenseVector<ValueType>& onePlusLtauP,lama::DenseVector<ValueType>& onePlusLtauS){
    
    lama::DenseVector<ValueType>& pWaveModulus=model.getPWaveModulus();
    lama::DenseVector<ValueType>& sWaveModulus=model.getSWaveModulus();
    
    lama::DenseVector<ValueType> temp(sWaveModulus.getDistributionPtr());
    lama::DenseVector<ValueType> temp2(sWaveModulus.getDistributionPtr());

    /* --------------------------------------- */
    /* Apply scaling for update of Sxx and Szz */
    /* --------------------------------------- */

    temp=2*sWaveModulus;
    temp.scale(onePlusLtauS);
    
    temp2=-1.0*pWaveModulus;
    temp2.scale(onePlusLtauP);
    
    temp += temp2; // = ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) );
    
    scaleStressHorizontalUpdate=pWaveModulus;
    scaleStressHorizontalUpdate.scale(onePlusLtauP); // = ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.invert();  // = 1 / ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.scale(temp); // = ( ( 2 * S-wave Modul * ( 1 + L * tauS) ) -  ( P-wave Modul * ( 1 + L * tauP) ) ) / ( ( P-wave Modul * ( 1 + L * tauP) )
    scaleStressHorizontalUpdate.scale(selectHorizontalUpdate); // set to zero everywhere besides the surface
    
    /* --------------------------------------- */
    /* Apply scaling for update of Rxx and Rzz */
    /* --------------------------------------- */

    temp=2*sWaveModulus;
    temp.scale(tauS);
    
    temp2=-1.0*pWaveModulus;
    temp2.scale(tauP);
    
    temp += temp2; // = ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP );
    
    scaleRelaxationHorizontalUpdate=pWaveModulus;
    scaleRelaxationHorizontalUpdate.scale(tauP); // = ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate.invert();  // = 1 / ( P-wave Modul * tauP )
    scaleRelaxationHorizontalUpdate.scale(temp); // = ( ( 2 * S-wave Modul * tauS ) -  ( P-wave Modul * tauP ) ) / ( ( P-wave Modul tauP) )
    scaleRelaxationHorizontalUpdate.scale(selectHorizontalUpdate); // set to zero everywhere besides the surface
    

}


/*! \brief Apply free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonatlDerivative Sum of horizontal velocity updates
 \param temp DenseVector which is used for temporary storage
 \param Sxx Sxx wavefield (stress)
 \param Syy Syy wavefield (stress)
 \param Szz Szz wavefield (stress)
 \param Rxx Rxx wavefield (relaxation)
 \param Ryy Ryy wavefield (relaxation)
 \param Rzz Rzz wavefield (relaxation)
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<ValueType>::apply(lama::Vector& sumHorizonatlDerivative, lama::Vector& temp, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz, lama::DenseVector<ValueType>& Rxx, lama::DenseVector<ValueType>& Ryy, lama::DenseVector<ValueType>& Rzz){
    
    
    /* Update the stress parameter at the free surface */
    temp=sumHorizonatlDerivative;
    temp.scale(scaleStressHorizontalUpdate);
    Sxx +=temp; // Apply horizontal update
    Szz +=temp; // Apply horizontal update
    
    Syy.scale(setSurfaceZero); // Set the free surface to zero

    
    /* Update relaxation parameter at the free surface */
    sumHorizonatlDerivative.scale(scaleRelaxationHorizontalUpdate);
    Rxx +=sumHorizonatlDerivative; // Apply horizontal update
    Rzz +=sumHorizonatlDerivative; // Apply horizontal update
    
    Ryy.scale(setSurfaceZero); // Set the free surface to zero

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
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dvisco<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH){
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Initialization of the free surface...\n" );
    
    active=true;
    
    Derivatives::FD3D<ValueType>* derivativeFreeSurface= static_cast<Derivatives::FD3D<ValueType>*>(&derivatives);
    derivativeFreeSurface->useFreeSurface=true;
    derivativeFreeSurface->calcDyfPressure(NX, NY, NZ, dist);
    derivativeFreeSurface->calcDyfVelocity(NX, NY, NZ, dist);
    derivativeFreeSurface->calcDybPressure(NX, NY, NZ, dist);
    derivativeFreeSurface->calcDybVelocity(NX, NY, NZ, dist);
    derivativeFreeSurface->DybPressure.scale(lama::Scalar(DT/DH));
    derivativeFreeSurface->DybVelocity.scale(lama::Scalar(DT/DH));
    derivativeFreeSurface->DyfPressure.scale(lama::Scalar(DT/DH));
    derivativeFreeSurface->DyfVelocity.scale(lama::Scalar(DT/DH));
    derivativeFreeSurface->Dyb.purge();
    derivativeFreeSurface->Dyf.purge();
    
    selectHorizontalUpdate.allocate(dist);
    selectHorizontalUpdate=0.0;
    
    setSurfaceZero.allocate(dist);
    setSurfaceZero=1.0;
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); /* get local indices based on used distribution */
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    
    /* Get write access to local part of scaleHorizontalUpdate */
    utilskernel::LArray<ValueType>* selectHorizontalUpdate_LA=&selectHorizontalUpdate.getLocalValues();
    hmemo::WriteAccess<ValueType> write_selectHorizontalUpdate(*selectHorizontalUpdate_LA);
    
    /* Get write access to local part of setSurfaceZero */
    utilskernel::LArray<ValueType>* setSurfaceZero_LA=&setSurfaceZero.getLocalValues();
    hmemo::WriteAccess<ValueType> write_setSurfaceZero(*setSurfaceZero_LA);
    
    KITGPI::Acquisition::Coordinates<ValueType> coordinateTransformation;
    
    IndexType rowGlobal;
    IndexType rowLocal;
    
    for(IndexType i=0; i<numLocalIndices; i++){
        
        rowGlobal=read_localIndices[i];
        rowLocal=dist->global2local(rowGlobal);
        
        /* Determine if the current grid point is located on the surface */
        if(coordinateTransformation.locatedOnSurface(rowGlobal,NX,NY,NZ)){
            
            /* Set horizontal update to +1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal]=1.0;
            
            /* Set vector at surface to zero  */
            write_setSurfaceZero[rowLocal]=0.0;
        }
        
    }
    read_localIndices.release();
    write_selectHorizontalUpdate.release();
    write_setSurfaceZero.release();
    HOST_PRINT( dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n" );
}
