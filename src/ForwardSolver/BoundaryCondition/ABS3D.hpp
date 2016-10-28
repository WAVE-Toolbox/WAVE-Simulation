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
			class ABS3D : public ABS<ValueType>
			{
			public:

				//! Default constructor
				ABS3D(){};

				//! Default destructor
				~ABS3D(){};

				void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm );

                void apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4);
                void apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4, lama::DenseVector<ValueType>& v5, lama::DenseVector<ValueType>& v6, lama::DenseVector<ValueType>& v7, lama::DenseVector<ValueType>& v8, lama::DenseVector<ValueType>& v9);
                
			private:
                
                lama::DenseVector<ValueType> damping; //!< Absorbing Coefficient vector

				void absorbing_coefficients( IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist );

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
 \param v4 DenseVector to apply damping boundary
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4){
    
    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
    v4.scale(damping);

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
 \param v6 DenseVector to apply damping boundary
 \param v7 DenseVector to apply damping boundary
 \param v8 DenseVector to apply damping boundary
 \param v9 DenseVector to apply damping boundary
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::apply(lama::DenseVector<ValueType>& v1, lama::DenseVector<ValueType>& v2, lama::DenseVector<ValueType>& v3, lama::DenseVector<ValueType>& v4, lama::DenseVector<ValueType>& v5, lama::DenseVector<ValueType>& v6, lama::DenseVector<ValueType>& v7, lama::DenseVector<ValueType>& v8, lama::DenseVector<ValueType>& v9){
    
    v1.scale(damping);
    v2.scale(damping);
    v3.scale(damping);
    v4.scale(damping);
    v5.scale(damping);
    v6.scale(damping);
    v7.scale(damping);
    v8.scale(damping);
    v9.scale(damping);

}

//! \brief Calculation of absorbing coefficients
/*!
 *
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param dist Distribution of the wavefield
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::absorbing_coefficients( IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist )
{
	IndexType FW=20;

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
	ValueType DAMPING=8.0;
	ValueType amp=0;
	ValueType coeff[FW];
	coeff[0]++;  //need a better solution "todo"
	ValueType a=0;

	amp=1.0-DAMPING/100.0;
	a=sqrt ( -log ( amp ) / ( ( FW ) * ( FW ) ) );

	for( IndexType j=0; j<FW; j++ ) {
		coeff[j]=exp ( - ( a*a* ( FW-j ) * ( FW-j ) ) );
	}
    
    Acquisition::Coordinates<ValueType> coordTransform;
    Acquisition::coordinate3D coordinate;

    IndexType coordinateMin;
    
	/* Set the values into the indice arrays and the value array */
	for( IndexType i=0; i<numLocalIndices; i++ ) {

		read_localIndices_temp=read_localIndices[i];

		coordinate=coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ );
		coordinate=coordTransform.edgeDistance(coordinate, NX, NY, NZ );

        coordinateMin=coordinate.min();
		if( coordinateMin < FW ) {
            write_damping[i]=coeff[coordinateMin];
		}

	}
    
    /* Release all read and write access */
	read_localIndices.release();
    write_damping.release();

}

//! \brief Initializsation of the absorbing coefficient matrix
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::ABS3D<ValueType>::init ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )
{

	HOST_PRINT ( comm, "Initialization of the Damping Boundary...\n" );

	absorbing_coefficients ( NX, NY, NZ, dist );
    
	damping.setContextPtr ( ctx );

	HOST_PRINT ( comm, "Finished with initialization of the Damping Boundary!\n\n" );

}






