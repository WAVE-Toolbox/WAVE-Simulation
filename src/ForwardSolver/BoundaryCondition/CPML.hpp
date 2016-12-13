#pragma once

namespace KITGPI {
    
    namespace ForwardSolver {
        
        namespace BoundaryCondition {
            
            //! \brief Abstract class for the calculation and application of cpml boundaries
            template<typename ValueType>
            class CPML
            {
            public:
                
                //! \brief Default constructor
                CPML(){};
                
                //! \brief Default destructor
                ~CPML(){};
		
                //! init CPML coefficient vectors and CPML memory variables
                virtual void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ,ValueType DT, IndexType DH, IndexType BoundaryWidth,bool useFreeSurface,Configuration::PMLVariables<ValueType> const &PMLVar)=0;
                    
	    protected:
		void resetVector ( lama::DenseVector<ValueType>& vector );

		void initVector ( lama::DenseVector<ValueType>& vector,hmemo::ContextPtr ctx, dmemo::DistributionPtr dist );
		    
		void SetCoeffCPML ( lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv,lama::DenseVector<ValueType>& a_half, lama::DenseVector<ValueType>& b_half,lama::DenseVector<ValueType>& kInv_half,IndexType coord,
				IndexType gdist, IndexType BoundaryWidth,Configuration::PMLVariables<ValueType> const&PMLVar,IndexType i, ValueType DT , ValueType DH );
		
		void ResetCoeffFreeSurface ( lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv,
						lama::DenseVector<ValueType>& a_half, lama::DenseVector<ValueType>& b_half,lama::DenseVector<ValueType>& kInv_half,
						IndexType i);
		
		/*inline*/ void applyCPML ( lama::Vector& Vec,lama::DenseVector<ValueType>& Psi,lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv );
		    
		    		lama::DenseVector<ValueType> psi_vxx;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vyx;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vzx;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vxy;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vyy;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vzy;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vxz;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vyz;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_vzz;//!< CPML memory Variable

				lama::DenseVector<ValueType> psi_sxx_x;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_sxy_x;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_sxz_x;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_sxy_y;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_syy_y;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_syz_y;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_sxz_z;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_syz_z;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_szz_z;//!< CPML memory Variable
				
				lama::DenseVector<ValueType> psi_p_x;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_p_y;//!< CPML memory Variable
				lama::DenseVector<ValueType> psi_p_z;//!< CPML memory Variable
				
				
				
				lama::DenseVector<ValueType> k_x;//!< CPML coefficient
				lama::DenseVector<ValueType> b_x;//!< CPML coefficient
				lama::DenseVector<ValueType> a_x;//!< CPML coefficient
				lama::DenseVector<ValueType> k_y;//!< CPML coefficient
				lama::DenseVector<ValueType> b_y;//!< CPML coefficient
				lama::DenseVector<ValueType> a_y;//!< CPML coefficient
				lama::DenseVector<ValueType> k_z;//!< CPML coefficient
				lama::DenseVector<ValueType> b_z;//!< CPML coefficient
				lama::DenseVector<ValueType> a_z;//!< CPML coefficient

				lama::DenseVector<ValueType> k_x_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> b_x_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> a_x_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> k_y_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> b_y_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> a_y_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> k_z_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> b_z_half;//!< CPML coefficient for staggered gridpoints
				lama::DenseVector<ValueType> a_z_half;//!< CPML coefficient for staggered gridpoints
				
				lama::DenseVector<ValueType> update_PmlTemp;//!< temporary vector for pml application
		
            };
        } /* end namespace BoundaryCondition  */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */

/*! \brief Reset a single Vector to zero.
*/
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::resetVector ( lama::DenseVector<ValueType>& vector )
{
	vector=0.0;
}

/*! \brief Intitialisation of a single vector.
*
* This method will set the context, allocate the the wavefield and set the field to zero.
*/
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::initVector ( lama::DenseVector<ValueType>& vector,hmemo::ContextPtr ctx, dmemo::DistributionPtr dist )
{
	vector.setContextPtr ( ctx );
	vector.allocate ( dist );

	resetVector ( vector );
}

/*! \brief set CPML coefficients
 * 
 * method to set cpml coefficients for a given gridpoint
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::SetCoeffCPML ( lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv,
										 lama::DenseVector<ValueType>& a_half, lama::DenseVector<ValueType>& b_half,lama::DenseVector<ValueType>& kInv_half,
										 IndexType coord,IndexType gdist, IndexType BoundaryWidth,Configuration::PMLVariables<ValueType> const &PMLVar, 
										 IndexType i, ValueType DT, ValueType DH )
{

	ValueType NPower=PMLVar.NPower;
	ValueType K_Max_Pml=PMLVar.KMaxCPML;
	
	ValueType RCoef=0.0008;
	
	ValueType alpha_max_Pml=2.0 * M_PI * ( PMLVar.CenterFrequencyCPML/2.0 );
	
	ValueType d0 = - ( NPower + 1 ) * PMLVar.VMaxCPML * log ( RCoef ) / ( 2.0 * BoundaryWidth*DH );
	
	
	ValueType PositionNorm=0.0;


	ValueType k_temp=0.0;
	ValueType b_temp=0.0;
	ValueType a_temp=0.0;

	ValueType alpha_prime=0.0;
	ValueType d=0.0;

	utilskernel::LArray<ValueType>* k_LA=&kInv.getLocalValues();
	hmemo::WriteAccess<ValueType> write_k ( *k_LA );

	utilskernel::LArray<ValueType>* b_LA=&b.getLocalValues();
	hmemo::WriteAccess<ValueType> write_b ( *b_LA );

	utilskernel::LArray<ValueType>* a_LA=&a.getLocalValues();
	hmemo::WriteAccess<ValueType> write_a ( *a_LA );

	utilskernel::LArray<ValueType>* k_half_LA=&kInv_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_k_half ( *k_half_LA );

	utilskernel::LArray<ValueType>* b_half_LA=&b_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_b_half ( *b_half_LA );

	utilskernel::LArray<ValueType>* a_half_LA=&a_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_a_half ( *a_half_LA );


	/* left boundary */
	if ( coord < BoundaryWidth ) {
		PositionNorm= ( ValueType ) ( BoundaryWidth-gdist ) /BoundaryWidth;
		d = d0 * pow ( PositionNorm,NPower );
		k_temp = 1.0 + ( K_Max_Pml - 1.0 ) * pow ( PositionNorm,NPower );
		alpha_prime = alpha_max_Pml * ( 1.0 - PositionNorm );
		b_temp = exp ( - ( d / k_temp + alpha_prime ) * DT );
		/* avoid division by zero outside the PML */
		if ( abs ( d ) > 1.0e-6 ) {
			a_temp = d * ( b_temp - 1.0 ) / ( k_temp * ( d + k_temp * alpha_prime ) );
		} else a_temp=0.0;
		write_a[i]=a_temp;
		write_b[i]= b_temp;
		write_k[i]=1/k_temp;


		/* right boundary
		 *starts with half point -> first point shiftet +1 to th right*/
	} else if ( gdist < BoundaryWidth-1 ) {
		PositionNorm= ( ValueType ) ( BoundaryWidth-gdist-1 ) /BoundaryWidth;
		d = d0 * pow ( PositionNorm,NPower );
		k_temp = 1.0 + ( K_Max_Pml - 1.0 ) * pow ( PositionNorm,NPower );
		alpha_prime = alpha_max_Pml * ( 1.0 - PositionNorm );
		b_temp = exp ( - ( d / k_temp + alpha_prime ) * DT );
		/* avoid division by zero outside the PML */
		if ( abs ( d ) > 1.0e-6 ) {
			a_temp = d * ( b_temp - 1.0 ) / ( k_temp * ( d + k_temp * alpha_prime ) );
		} else a_temp=0.0;
		write_a[i]=a_temp;
		write_b[i]= b_temp;
		write_k[i]=1/k_temp;

	}



	/* half points */
	PositionNorm= ( ValueType ) ( BoundaryWidth-gdist-0.5 ) /BoundaryWidth;
	d = d0 * pow ( PositionNorm,NPower );
	k_temp = 1.0 + ( K_Max_Pml - 1.0 ) * pow ( PositionNorm,NPower );
	alpha_prime = alpha_max_Pml * ( 1.0 - PositionNorm );
	b_temp = exp ( - ( d / k_temp + alpha_prime ) * DT );
	/* avoid division by zero outside the PML */
	if ( abs ( d ) > 1.0e-6 ) {
		a_temp = d * ( b_temp - 1.0 ) / ( k_temp * ( d + k_temp * alpha_prime ) );
	} else a_temp=0.0;
	write_a_half[i]= a_temp;
	write_b_half[i]= b_temp;
	write_k_half[i]=1/k_temp;

	write_k.release();
	write_a.release();
	write_b.release();
	write_k_half.release();
	write_a_half.release();
	write_b_half.release();


}

/*! \brief Reset CPML coefficients for free surface
 * 
 * method to set cpml coefficients for a given gridpoint
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::ResetCoeffFreeSurface ( lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv,
										 lama::DenseVector<ValueType>& a_half, lama::DenseVector<ValueType>& b_half,lama::DenseVector<ValueType>& kInv_half,
										 IndexType i)
{


	utilskernel::LArray<ValueType>* k_LA=&kInv.getLocalValues();
	hmemo::WriteAccess<ValueType> write_k ( *k_LA );

	utilskernel::LArray<ValueType>* b_LA=&b.getLocalValues();
	hmemo::WriteAccess<ValueType> write_b ( *b_LA );

	utilskernel::LArray<ValueType>* a_LA=&a.getLocalValues();
	hmemo::WriteAccess<ValueType> write_a ( *a_LA );

	utilskernel::LArray<ValueType>* k_half_LA=&kInv_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_k_half ( *k_half_LA );

	utilskernel::LArray<ValueType>* b_half_LA=&b_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_b_half ( *b_half_LA );

	utilskernel::LArray<ValueType>* a_half_LA=&a_half.getLocalValues();
	hmemo::WriteAccess<ValueType> write_a_half ( *a_half_LA );

		write_a[i] = 0.0;
		write_b[i] = 0.0;
		write_k[i] = 1.0;

		write_a_half[i] = 0.0;
		write_b_half[i] = 0.0;
		write_k_half[i] = 1.0;

	write_k.release();
	write_a.release();
	write_b.release();
	write_k_half.release();
	write_a_half.release();
	write_b_half.release();
}




/*! \brief Application of the CPML
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param Vec DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 \param a CPML coefficient a
 \param b CPML coefficient b
 \param kInv reciprocal of CPML coefficient k
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML<ValueType>::applyCPML ( lama::Vector& Vec,lama::DenseVector<ValueType>& Psi,lama::DenseVector<ValueType>& a, lama::DenseVector<ValueType>& b,lama::DenseVector<ValueType>& kInv )
{

	update_PmlTemp=Vec;
	Psi.scale ( b );
	Psi+=update_PmlTemp.scale ( a );
	Vec.scale ( kInv );
	Vec+= Psi;
}