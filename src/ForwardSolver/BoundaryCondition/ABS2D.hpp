#pragma once

#include "ABS.hpp"
#include "../../Common/HostPrint.hpp"

namespace KITGPI
{
    
    namespace ForwardSolver
    {
        
        namespace BoundaryCondition
        {
            
            //! \brief Class for the calculation of the Absorbing Coefficient matrix for 3-D FD Simulations
            /*!
             * Calculation of the absorbing coefficient matrix for an equidistand grid
             *
             */
            template<typename ValueType>
            class ABS2D : public ABS<ValueType>
            {
            public:
                
                //! Default constructor
                ABS2D(){};
                
                //! Default destructor
                ~ABS2D(){};
                
                void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface);
                
                void apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3);
                void apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4, lama::DenseVector<ValueType>& v5);
                
            private:
                
                lama::DenseVector<ValueType> damping; //!< Absorbing Coefficient vector
                
            };
        } /* end namespace BoundaryCondition */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Application of the damping boundary
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param v1 DenseVector to apply damping boundary
 \param v2 DenseVector to apply damping boundary
 \param v3 DenseVector to apply damping boundary
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3){
    
    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
    
}

/*! \brief Application of the damping boundary
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param v1 DenseVector to apply damping boundary
 \param v2 DenseVector to apply damping boundary
 \param v3 DenseVector to apply damping boundary
 \param v4 DenseVector to apply damping boundary
 \param v5 DenseVector to apply damping boundary
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4, lama::DenseVector<ValueType>& v5){
    
    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
    v4.scale(damping);
    v5.scale(damping);
}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param BoundaryWidth Width of damping boundary
 \param DampingCoeff Damping coefficient
 \param useFreeSurface Bool if free surface is in use
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS2D<ValueType>::init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ,IndexType BoundaryWidth, ValueType DampingCoeff, bool useFreeSurface)
{
    
    HOST_PRINT ( dist->getCommunicatorPtr(), "Initialization of the Damping Boundary...\n" );
    
    dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);
    
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp; // Temporary storage, so we do not have to access the array
    
    /* Distributed vectors */
    damping.allocate(dist); // Vector to set elements on surface to zero
    damping=1.0;
    
    /* Get write access to local part of setSurfaceZero */
    utilskernel::LArray<ValueType>* damping_LA=&damping.getLocalValues();
    hmemo::WriteAccess<ValueType> write_damping(*damping_LA);
    
    // calculate damping function
    ValueType amp=0;
    ValueType coeff[BoundaryWidth];
    ValueType a=0;
    
    amp=1.0-DampingCoeff/100.0;
    a=sqrt ( -log ( amp ) / ( ( BoundaryWidth ) * ( BoundaryWidth ) ) );
    
    for( IndexType j=0; j<BoundaryWidth; j++ ) {
        coeff[j]=exp ( - ( a*a* ( BoundaryWidth-j ) * ( BoundaryWidth-j ) ) );
    }
    
    Acquisition::Coordinates<ValueType> coordTransform;
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D coordinatedist;
    
    IndexType coordinateMin=0;
    IndexType temp;
    
    /* Set the values into the indice arrays and the value array */
    for( IndexType i=0; i<numLocalIndices; i++ ) {
        
        read_localIndices_temp=read_localIndices[i];
        
        coordinate=coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ );
        coordinatedist=coordTransform.edgeDistance(coordinate, NX, NY, NZ );
        
        temp=0;
        if(coordinatedist.x<coordinatedist.y){
            temp=coordinatedist.x;
        } else {
            temp=coordinatedist.y;
        }
        coordinateMin=temp;
        if( coordinateMin < BoundaryWidth ) {
            write_damping[i]=coeff[coordinateMin];
        }
        
        if(useFreeSurface){
            if (coordinate.y < BoundaryWidth) {
                write_damping[i]=1.0;
                
                if ((coordinatedist.x < BoundaryWidth)){
                    write_damping[i]=coeff[coordinatedist.x];
                }
                
            }
        }
        
    }
    
    /* Release all read and write access */
    read_localIndices.release();
    write_damping.release();
    
    damping.setContextPtr ( ctx );
    
    HOST_PRINT ( dist->getCommunicatorPtr(), "Finished with initialization of the Damping Boundary!\n\n" );
    
}
