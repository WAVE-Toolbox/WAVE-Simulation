#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../Derivatives/FDTD3D.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurface.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D elastic free surface
            template<typename ValueType>
            class FreeSurfaceElastic : public FreeSurface<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurfaceElastic(){};
                
                virtual ~FreeSurfaceElastic()=0;
                
                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH);
                
                void setModelparameter(Modelparameter::Modelparameter<ValueType>& model);
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz);
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy);

            protected:
                
                using FreeSurface<ValueType>::active;
                
                lama::DenseVector<ValueType> setSurfaceZero; //!< Vector, which sets the wavefields at the surface to zero
                lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                lama::DenseVector<ValueType> scaleHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */

//! Default destructor
template<typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::~FreeSurfaceElastic(){};

/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType>& model){
    
    lama::DenseVector<ValueType>& pWaveModulus=model.getPWaveModulus();
    lama::DenseVector<ValueType>& sWaveModulus=model.getSWaveModulus();
    
    lama::DenseVector<ValueType> temp(sWaveModulus.getDistributionPtr());
    
    temp=2*sWaveModulus-pWaveModulus;
    
    scaleHorizontalUpdate=pWaveModulus;
    scaleHorizontalUpdate.invert();
    scaleHorizontalUpdate.scale(temp);
    scaleHorizontalUpdate.scale(selectHorizontalUpdate);

}


/*! \brief Apply free surface condition during time stepping for 2D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy){
    
    SCAI_ASSERT_DEBUG( active , " FreeSurface is not active " );
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    
    Sxx +=sumHorizonalDerivative;
    
    Syy.scale(setSurfaceZero);
    
}

/*! \brief Apply free surface condition during time stepping for 3D simulations
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sumHorizonalDerivative Sum of horizontal velocity updates
 \param Sxx Sxx wavefield
 \param Syy Syy wavefield
 \param Szz Szz wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz){
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    
    Sxx +=sumHorizonalDerivative;
    Szz +=sumHorizonalDerivative;
    
    Syy.scale(setSurfaceZero);
    
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH){
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Initialization of the free surface...\n" );
    
    active=true;
    
    derivatives.useFreeSurface=true;
    derivatives.calcDyfPressure(NX, NY, NZ, dist);
    derivatives.calcDyfVelocity(NX, NY, NZ, dist);
    derivatives.calcDybPressure(NX, NY, NZ, dist);
    derivatives.calcDybVelocity(NX, NY, NZ, dist);
    derivatives.DybPressure.scale(lama::Scalar(DT/DH));
    derivatives.DybVelocity.scale(lama::Scalar(DT/DH));
    derivatives.DyfPressure.scale(lama::Scalar(DT/DH));
    derivatives.DyfVelocity.scale(lama::Scalar(DT/DH));
    derivatives.Dyb.purge();
    derivatives.Dyf.purge();
    
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
            
            /* Set horizontal update to 1 at the surface and leave it zero else */
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
