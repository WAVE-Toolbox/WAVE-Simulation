#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurface.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D elastic free surface
            template<typename ValueType>
            class FreeSurface3Delastic : public FreeSurface<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface3Delastic(){};
                
                //! Default destructor
                ~FreeSurface3Delastic(){};
                
                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ);
                
                void setModelparameter(Modelparameter::Modelparameter<ValueType>& model);
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz);
                
            private:
                
                using FreeSurface<ValueType>::active;
                
                lama::DenseVector<ValueType> setSurfaceZero; //!< Vector, which sets the wavefields at the surface to zero
                lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                lama::DenseVector<ValueType> scaleHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */



/*! \brief Scale horizontal update with model parameter
 *
 *
 \param model which is used during forward modelling
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType>& model){
    
    lama::DenseVector<ValueType>& lambda=model.getLambda();
    lama::DenseVector<ValueType>& mu=model.getMu();
    
    scaleHorizontalUpdate=lambda+2*mu;
    scaleHorizontalUpdate.invert();
    scaleHorizontalUpdate.scale(lambda);
    scaleHorizontalUpdate.scale(selectHorizontalUpdate);

}




/*! \brief Apply free surface condition during time stepping
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
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& /*Syy*/, lama::DenseVector<ValueType>& Szz){
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    
    Sxx +=sumHorizonalDerivative;
    Szz +=sumHorizonalDerivative;
    
}

/*! \brief Initialitation of the free surface
 *
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param NX Number of grid points in X-Direction
 \param NY Number of grid points in Y-Direction (Depth)
 \param NZ Number of grid points in Z-Direction
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& /*derivatives*/, IndexType NX, IndexType NY, IndexType NZ){
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Initialization of the free surface...\n" );
    
    active=true;
    
    selectHorizontalUpdate.allocate(dist);
    selectHorizontalUpdate=0.0;
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); /* get local indices based on used distribution */
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    
    /* Get write access to local part of scaleHorizontalUpdate */
    utilskernel::LArray<ValueType>* selectHorizontalUpdate_LA=&selectHorizontalUpdate.getLocalValues();
    hmemo::WriteAccess<ValueType> write_selectHorizontalUpdate(*selectHorizontalUpdate_LA);
    
    KITGPI::Acquisition::Coordinates<ValueType> coordinateTransformation;
    
    IndexType rowGlobal;
    IndexType rowLocal;
    
    for(IndexType i=0; i<numLocalIndices; i++){
        
        rowGlobal=read_localIndices[i];
        rowLocal=dist->global2local(rowGlobal);
        
        /* Determine if the current grid point is located on the surface */
        if(coordinateTransformation.locatedOnSurface(rowGlobal,NX,NY,NZ)){
            
            /* Set horizontal update to -1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal]=-1.0;
            
        }
        
    }
    read_localIndices.release();
    write_selectHorizontalUpdate.release();

    HOST_PRINT( dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n" );
}
