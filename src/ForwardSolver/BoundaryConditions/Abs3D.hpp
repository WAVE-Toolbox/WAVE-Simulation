#pragma once

#include "Abs.hpp"
#include "../../Common/HostPrint.hpp"

namespace KITGPI
{

	namespace ForwardSolver
	{

		//! \brief Abs namespace
		namespace Abs
		{


			//! \brief Class for the calculation of the Absorbing Coefficient matrix for 3-D FD Simulations
			/*!
			 * Calculation of the absorbing coefficient matrix for an equidistand grid
			 *
			 */
			template<typename ValueType>
			class Abs3D : public Abs<ValueType>
			{
			public:

				//! Default constructor
				Abs3D() {};

				//! Default destructor
				~Abs3D() {};

				Abs3D ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm );
				Abs3D ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm );

				void initializeMatrix ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm );

				lama::CSRSparseMatrix<ValueType>& getAbsCoeff();

			private:

				lama::CSRSparseMatrix<ValueType> AbsCoeff; //!< Absorbing Coefficient Matrix

				void absorbing_coefficients ( IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist );

			};
		} /* end namespace Abs */
	} /* end namespace ForwardSolver */
} /* end namespace KITGPI */


//! \brief Constructor to support Configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template<typename ValueType>
KITGPI::ForwardSolver::Abs::Abs3D<ValueType>::Abs3D ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx, Configuration::Configuration<ValueType> config, dmemo::CommunicatorPtr comm )
{
	Abs<ValueType>::initializeMatrix ( dist,ctx, config, comm );
}

//! \brief Constructor of the absorbing coefficient matrix
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
KITGPI::ForwardSolver::Abs::Abs3D<ValueType>::Abs3D ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )
{
	initializeMatrix ( dist, ctx, NX, NY,NZ, comm );
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
void KITGPI::ForwardSolver::Abs::Abs3D<ValueType>::absorbing_coefficients ( IndexType NX, IndexType NY, IndexType NZ, dmemo::DistributionPtr dist )
{
	IndexType FW=20;

	dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();

	/* Get local "global" indices */
	hmemo::HArray<IndexType> localIndices;
	dist->getOwnedIndexes ( localIndices );

	IndexType N=NX*NY*NZ;

	IndexType numLocalIndices=localIndices.size(); //< Number of local indices

	/* Add the number of non-diagonal values */
	hmemo::ReadAccess<IndexType> read_localIndices ( localIndices ); //< Get read access to localIndices
	IndexType read_localIndices_temp; //< Temporary storage, so we do not have to access the array

	/* Allocate local part to create local CSR storage*/
	hmemo::HArray<ValueType> AbsCoeff_valuesLocal ( numLocalIndices );
	hmemo::HArray<IndexType> AbsCoeff_csrJALocal ( numLocalIndices );
	hmemo::HArray<IndexType> AbsCoeff_csrIALocal ( numLocalIndices+1 );


	/* Get WriteAccess to local part */
	hmemo::WriteAccess<IndexType> AbsCoeff_write_csrJALocal ( AbsCoeff_csrJALocal );
	hmemo::WriteAccess<IndexType> AbsCoeff_write_csrIALocal ( AbsCoeff_csrIALocal );
	hmemo::WriteAccess<ValueType> AbsCoeff_write_valuesLocal ( AbsCoeff_valuesLocal );


	/* Set some counters to create the CSR Storage */
	IndexType AbsCoeff_countJA=0;
	IndexType AbsCoeff_countIA=0;
	AbsCoeff_write_csrIALocal[0]=0;
	AbsCoeff_countIA++;

	// calculate damping function

	ValueType DAMPING=8.0;
	ValueType amp=0;
	ValueType coeff[FW];
	coeff[0]++;  //need a better solution "todo"
	ValueType a=0;

	amp=1.0-DAMPING/100.0;
	a=sqrt ( -log ( amp ) / ( ( FW ) * ( FW ) ) );

	for ( IndexType j=0; j<FW; j++ ) {
		coeff[j]=exp ( - ( a*a* ( FW-j ) * ( FW-j ) ) );
	}


	/* Set the values into the indice arrays and the value array */
	for ( IndexType i=0; i<numLocalIndices; i++ ) {

		read_localIndices_temp=read_localIndices[i];
		//      read_localIndices_temp_plusOne=read_localIndices_temp+1;

		Acquisition::coordinate3D coord;
		Acquisition::distance3D distance;
		Acquisition::Coordinates<ValueType> coordclass;
		coord=coordclass.index2coordinate ( read_localIndices_temp, NX, NY, NZ );
		distance=coordclass.edgeDistance ( coord, NX, NY, NZ );



		if ( distance.min < FW ) {
			AbsCoeff_write_csrJALocal[AbsCoeff_countJA]=read_localIndices_temp;
			AbsCoeff_write_valuesLocal[AbsCoeff_countJA]=coeff[distance.min];
			AbsCoeff_countJA++;
		} else {

			AbsCoeff_write_csrJALocal[AbsCoeff_countJA]=read_localIndices_temp;
			AbsCoeff_write_valuesLocal[AbsCoeff_countJA]=1.0;
			AbsCoeff_countJA++;
		}


		AbsCoeff_write_csrIALocal[AbsCoeff_countIA]=AbsCoeff_countJA;
		AbsCoeff_countIA++;


	}
//	std::cout << "JA " << AbsCoeff_countJA << "  IA  " <<AbsCoeff_countIA << std::endl;

//     /* Release all read and write access */
	read_localIndices.release();

	AbsCoeff_write_csrJALocal.release();
	AbsCoeff_write_csrIALocal.release();
	AbsCoeff_write_valuesLocal.release();

	/* Create local CSR storage of Matrix A, than create distributed CSR matrix A */
	lama::CSRStorage<ValueType> AbsCoeff_LocalCSR ( numLocalIndices,N,numLocalIndices,AbsCoeff_csrIALocal,AbsCoeff_csrJALocal,AbsCoeff_valuesLocal );
	AbsCoeff.assign ( AbsCoeff_LocalCSR,dist,dist );
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
void KITGPI::ForwardSolver::Abs::Abs3D<ValueType>::initializeMatrix ( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ, dmemo::CommunicatorPtr comm )
{


	HOST_PRINT ( comm, "Initialization of the matrix AbsCoeff\n" );

	absorbing_coefficients ( NX, NY, NZ, dist );
	HOST_PRINT ( comm, "Matrix AbsCoeff finished.\n" );

	AbsCoeff.setContextPtr ( ctx );



	HOST_PRINT ( comm, "Finished with initialization of the matrix!\n" );


}

//! \brief Getter method for derivative matrix A
template<typename ValueType>
lama::CSRSparseMatrix<ValueType>& KITGPI::ForwardSolver::Abs::Abs3D<ValueType>::getAbsCoeff()
{
	return ( AbsCoeff );
}






