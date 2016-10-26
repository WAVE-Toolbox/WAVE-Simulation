#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            template<typename ValueType>
            class FreeSurface3Delastic
            {
            public:
                
                //! Default constructor
                FreeSurface3Delastic():active(false){};
                
                //! Default destructor
                ~FreeSurface3Delastic(){};
                
                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ);
                
                void setModelparameter(Modelparameter::Modelparameter<ValueType>& model);
                
                lama::CSRSparseMatrix<ValueType>& getDybPressure();
                
                bool getActive();
                
                void apply(lama::Vector& sumHorizonatlDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz);
                
            private:
                
                bool active;
                
                lama::DenseVector<ValueType> setSurfaceZero;
                lama::DenseVector<ValueType> selectHorizontalUpdate;

                lama::DenseVector<ValueType> scaleHorizontalUpdate;
                
                lama::CSRSparseMatrix<ValueType> DybPressure;
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::setModelparameter(Modelparameter::Modelparameter<ValueType>& model){
    
    lama::DenseVector<ValueType>& lambda=model.getLambda();
    lama::DenseVector<ValueType>& mu=model.getMu();
    
    scaleHorizontalUpdate=lambda+2*mu;
    scaleHorizontalUpdate.invert();
    scaleHorizontalUpdate.scale(lambda);
    scaleHorizontalUpdate.scale(selectHorizontalUpdate);

}

template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::getDybPressure(){
    if(!active){
        COMMON_THROWEXCEPTION(" Free Surface is not active ! ")
    }
    return(DybPressure);
}


template<typename ValueType>
bool KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::getActive(){
    return(active);
}

template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::apply(lama::Vector& sumHorizonalDerivative, lama::DenseVector<ValueType>& Sxx, lama::DenseVector<ValueType>& Syy, lama::DenseVector<ValueType>& Szz){
    
    /* Apply horizontal update, which replaces the vertical one */
    sumHorizonalDerivative.scale(scaleHorizontalUpdate);
    Sxx +=sumHorizonalDerivative;
    Szz +=sumHorizonalDerivative;
    
    /* Set the elements on the surface to zero */
    Syy.scale(setSurfaceZero);
    
}

template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Delastic<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ){
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Initialization of the free surface...\n" );
    
    active=true;
    
    lama::CSRSparseMatrix<ValueType>& DybVelocity=derivatives.getE();
    DybPressure=DybVelocity;
    
    IndexType size_vec=dist->getGlobalSize();
    
    /* Local vectors */
    lama::DenseVector<ValueType> zeroRowLocal(size_vec,0.0); // Zero vector
    lama::DenseVector<ValueType> modfiyRowLocal(size_vec,0.0); // Vector to manipulate row content of derivative matrix
    
    /* Distributed vectors */
    setSurfaceZero.allocate(dist);
    setSurfaceZero=1.0; // Vector to set elements on surface to zero
    
    selectHorizontalUpdate.allocate(dist);
    selectHorizontalUpdate=0.0;
    
    hmemo::HArray<ValueType> zeroRowHArray(size_vec,0.0);
    hmemo::HArray<ValueType> modfiyRowHArray(size_vec,0.0);
    
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
            
            /* Set vertical updates of the strain to zero at the surface */
            DybPressure.setLocalRow(zeroRowHArray,rowLocal,scai::utilskernel::binary::BinaryOp::MULT);
            
            /* Set horizontal update to -1 at the surface and leave it zero else */
            write_selectHorizontalUpdate[rowLocal]=-1.0;
            
            /* Set elements at the surface to zero */
            write_setSurfaceZero[rowLocal]=0.0;
            
            /* Modify vertical derivative matrix to apply imaging condition */
            modfiyRowHArray.init(size_vec,0.0);
            hmemo::WriteAccess<ValueType> write_modfiyRowHArray(modfiyRowHArray);
            write_modfiyRowHArray[rowGlobal]=2.0;
            write_modfiyRowHArray.release();
            DybVelocity.setLocalRow(modfiyRowHArray,rowLocal,scai::utilskernel::binary::BinaryOp::MULT);
            
        }
        
    }
    read_localIndices.release();
    write_selectHorizontalUpdate.release();
    write_setSurfaceZero.release();
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n" );
}
