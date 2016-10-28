#pragma once

#include "../Derivatives/Derivatives.hpp"
#include "../../Common/HostPrint.hpp"
#include "FreeSurface.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition {
            
            //! \brief 3-D acoustic free surface
            template<typename ValueType>
            class FreeSurface3Dacoustic : public FreeSurface<ValueType>
            {
            public:
                
                //! Default constructor
                FreeSurface3Dacoustic(){};
                
                //! Default destructor
                ~FreeSurface3Dacoustic(){};
                
                void init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ);
                
                void apply(lama::DenseVector<ValueType>& p);
                
            private:
                
                using FreeSurface<ValueType>::active;
                
                lama::DenseVector<ValueType> setSurfaceZero; //!< Vector, which sets the wavefields at the surface to zero
                
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */



/*! \brief Apply free surface condition during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param p p wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dacoustic<ValueType>::apply(lama::DenseVector<ValueType>& p){
    
    /* Set the elements on the surface to zero */
    p.scale(setSurfaceZero);
    
}

/*! \brief Initialitation of the free surface
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param NX Number of grid points in X-Direction
 \param NY Number of grid points in Y-Direction (Depth)
 \param NZ Number of grid points in Z-Direction
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurface3Dacoustic<ValueType>::init(dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType>& derivatives, IndexType NX, IndexType NY, IndexType NZ){
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Initialization of the free surface...\n" );
    
    active=true;
    
    lama::CSRSparseMatrix<ValueType>& Dyb=derivatives.getE();
    
    IndexType size_vec=dist->getGlobalSize();
    
    /* Local vectors */
    lama::DenseVector<ValueType> zeroRowLocal(size_vec,0.0); // Zero vector
    //lama::DenseVector<ValueType> modfiyRowLocal(size_vec,0.0); // Vector to manipulate row content of derivative matrix
    
    /* Distributed vectors */
    setSurfaceZero.allocate(dist); // Vector to set elements on surface to zero
    setSurfaceZero=1.0;
    
    hmemo::HArray<ValueType> zeroRowHArray(size_vec,0.0);
    //hmemo::HArray<ValueType> modfiyRowHArray(size_vec,0.0);
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); /* get local indices based on used distribution */
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    
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
            Dyb.setLocalRow(zeroRowHArray,rowLocal,scai::utilskernel::binary::BinaryOp::MULT);
            
            /* Set elements at the surface to zero */
            write_setSurfaceZero[rowLocal]=0.0;
            
            /* Modify vertical derivative matrix to apply imaging condition */
            //modfiyRowHArray.init(size_vec,0.0);
            //hmemo::WriteAccess<ValueType> write_modfiyRowHArray(modfiyRowHArray);
            //write_modfiyRowHArray[rowGlobal]=2.0;
            //write_modfiyRowHArray.release();
            //DybVelocity.setLocalRow(modfiyRowHArray,rowLocal,scai::utilskernel::binary::BinaryOp::MULT);
            
        }
        
    }
    read_localIndices.release();
    write_setSurfaceZero.release();
    
    HOST_PRINT( dist->getCommunicatorPtr(), "Finished initializing of the free surface\n\n" );
}
