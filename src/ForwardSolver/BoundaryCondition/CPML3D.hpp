#pragma once

#include "CPML.hpp"
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
			class CPML3D : public CPML<ValueType>
			{
			public:

				//! Default constructor
				CPML3D(){};

				//! Default destructor
				~CPML3D(){};

				void init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ,ValueType DT,IndexType DH, IndexType BoundaryWidth, bool useFreeSurface,Configuration::PMLVariables<ValueType> &PMLVar);
				
				
				
			void applyHalfX(lama::Vector& VecX,lama::DenseVector<ValueType>& Psi);
			void applyX(lama::Vector& VecX,lama::DenseVector<ValueType>& Psi);
			void applyHalfY(lama::Vector& VecY,lama::DenseVector<ValueType>& Psi);
			void applyY(lama::Vector& VecY,lama::DenseVector<ValueType>& Psi);
			void applyHalfZ(lama::Vector& VecZ,lama::DenseVector<ValueType>& Psi);
			void applyZ(lama::Vector& VecZ,lama::DenseVector<ValueType>& Psi);
			
			private:
				lama::DenseVector<ValueType> k_x;
				lama::DenseVector<ValueType> b_x;
				lama::DenseVector<ValueType> a_x;
				lama::DenseVector<ValueType> k_y;
				lama::DenseVector<ValueType> b_y;
				lama::DenseVector<ValueType> a_y;
				lama::DenseVector<ValueType> k_z;
				lama::DenseVector<ValueType> b_z;
				lama::DenseVector<ValueType> a_z;
				
				lama::DenseVector<ValueType> k_x_half;
				lama::DenseVector<ValueType> b_x_half;
				lama::DenseVector<ValueType> a_x_half;
				lama::DenseVector<ValueType> k_y_half;
				lama::DenseVector<ValueType> b_y_half;
				lama::DenseVector<ValueType> a_y_half;
				lama::DenseVector<ValueType> k_z_half;
				lama::DenseVector<ValueType> b_z_half;
				lama::DenseVector<ValueType> a_z_half;
				
				lama::DenseVector<ValueType> update_PmlTemp;

			};
		} /* end namespace BoundaryCondition */
	} /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Application of the CPML in x direction on half grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyHalfX(lama::Vector& VecX,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecX;
	Psi=Psi.scale(b_x_half)+update_PmlTemp.scale(a_x_half);
	VecX = VecX.scale(k_x_half) + Psi;

}

/*! \brief Application of the CPML in x direction on full grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyX(lama::Vector& VecX,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecX;
	Psi=Psi.scale(b_x)+update_PmlTemp.scale(a_x);
	VecX = VecX.scale(k_x) + Psi;

}

/*! \brief Application of the CPML in y direction on half grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyHalfY(lama::Vector& VecY,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecY;
	Psi=Psi.scale(b_y_half)+update_PmlTemp.scale(a_y_half);
	VecY = VecY.scale(k_y_half) + Psi;

}

/*! \brief Application of the CPML in y direction on full grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyY(lama::Vector& VecY,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecY;
	Psi=Psi.scale(b_y)+update_PmlTemp.scale(a_y);
	VecY = VecY.scale(k_y) + Psi;

}

/*! \brief Application of the CPML in z direction on half grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyHalfZ(lama::Vector& VecZ,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecZ;
	Psi=Psi.scale(b_z_half)+update_PmlTemp.scale(a_z_half);
	VecZ = VecZ.scale(k_z_half) + Psi;

}

/*! \brief Application of the CPML in z direction on full grid points
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param update DenseVector (derivate component) to apply pml
 \param Psi CPML Wavefield Vector
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::applyZ(lama::Vector& VecZ,lama::DenseVector<ValueType>& Psi){
    
	update_PmlTemp=VecZ;
	Psi=Psi.scale(b_z)+update_PmlTemp.scale(a_z);
	VecZ = VecZ.scale(k_z) + Psi;

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
 \param useFreeSurface Bool if free surface is in use
 */
template<typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::CPML3D<ValueType>::init(dmemo::DistributionPtr dist, hmemo::ContextPtr ctx,IndexType NX, IndexType NY, IndexType NZ,ValueType DT, IndexType DH, IndexType BoundaryWidth, bool useFreeSurface,Configuration::PMLVariables<ValueType> &PMLVar)
{
	
	HOST_PRINT ( dist->getCommunicatorPtr(), "Initialization of the PMl Coefficients...\n" );
    
    if(useFreeSurface){
        COMMON_THROWEXCEPTION(" Free Surface and CPML boundary are not implemented for simultaneous usage ! ")
    }
    
    ValueType NPower=PMLVar.NPower;
    ValueType K_Max_Pml=PMLVar.KMaxCPML;
    
    ValueType RCoef=0.0008;
    
    ValueType alpha_max_Pml=2.0 * M_PI * (PMLVar.CenterFrequencyCPML/2.0);
    
    ValueType d0;
    ValueType PositionNorm;
    
   d0 = - (NPower + 1) * PMLVar.VMaxCPML * log(RCoef) / (2.0 * BoundaryWidth*DH);
    
    ValueType k_temp=0.0;
    ValueType b_temp=0.0;
    ValueType a_temp=0.0;
   
    ValueType alpha_prime=0.0;
    ValueType d=0.0;
    
    dmemo::CommunicatorPtr comm=dist->getCommunicatorPtr();
    
    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);
    
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp; // Temporary storage, so we do not have to access the array
    
    update_PmlTemp.allocate(dist); 
    update_PmlTemp=0;
    
    
    /* Distributed vectors */
    k_x.allocate(dist); 
    k_y.allocate(dist); 
    k_z.allocate(dist); 
    b_x.allocate(dist); 
    b_y.allocate(dist); 
    b_z.allocate(dist);
    a_x.allocate(dist);
    a_y.allocate(dist);
    a_z.allocate(dist);
    
    k_x_half.allocate(dist); 
    k_y_half.allocate(dist); 
    k_z_half.allocate(dist); 
    b_x_half.allocate(dist); 
    b_y_half.allocate(dist); 
    b_z_half.allocate(dist);
    a_x_half.allocate(dist);
    a_y_half.allocate(dist);
    a_z_half.allocate(dist);
    
    k_x=1.0;
    k_y=1.0;
    k_z=1.0;
    b_x=0.0; 
    b_y=0.0;
    b_z=0.0;
    a_x=0.0;
    a_y=0.0;
    a_z=0.0;
    
    k_x_half=1.0;
    k_y_half=1.0;
    k_z_half=1.0;
    b_x_half=0.0;
    b_y_half=0.0;
    b_z_half=0.0;
    a_x_half=0.0;
    a_y_half=0.0;
    a_z_half=0.0;
    
    
    /* Get write access to local part of setSurfaceZero */
  utilskernel::LArray<ValueType>* k_x_LA=&k_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_x(*k_x_LA);
   
  utilskernel::LArray<ValueType>* b_x_LA=&b_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_x(*b_x_LA);

  utilskernel::LArray<ValueType>* a_x_LA=&a_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_x(*a_x_LA);

  utilskernel::LArray<ValueType>* k_y_LA=&k_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_y(*k_y_LA);

  utilskernel::LArray<ValueType>* b_y_LA=&b_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_y(*b_y_LA);

  utilskernel::LArray<ValueType>* a_y_LA=&a_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_y(*a_y_LA);

  utilskernel::LArray<ValueType>* k_z_LA=&k_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_z(*k_z_LA);

  utilskernel::LArray<ValueType>* b_z_LA=&b_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_z(*b_z_LA);

  utilskernel::LArray<ValueType>* a_z_LA=&a_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_z(*a_z_LA);
  
  
  utilskernel::LArray<ValueType>* k_x_half_LA=&k_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_x_half(*k_x_half_LA);

  utilskernel::LArray<ValueType>* b_x_half_LA=&b_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_x_half(*b_x_half_LA);

  utilskernel::LArray<ValueType>* a_x_half_LA=&a_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_x_half(*a_x_half_LA);

  utilskernel::LArray<ValueType>* k_y_half_LA=&k_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_y_half(*k_y_half_LA);

  utilskernel::LArray<ValueType>* b_y_half_LA=&b_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_y_half(*b_y_half_LA);

  utilskernel::LArray<ValueType>* a_y_half_LA=&a_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_y_half(*a_y_half_LA);

  utilskernel::LArray<ValueType>* k_z_half_LA=&k_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_z_half(*k_z_half_LA);

  utilskernel::LArray<ValueType>* b_z_half_LA=&b_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_z_half(*b_z_half_LA);

  utilskernel::LArray<ValueType>* a_z_half_LA=&a_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_z_half(*a_z_half_LA);

  
    Acquisition::Coordinates<ValueType> coordTransform;
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;

      for( IndexType i=0; i<numLocalIndices; i++ ) {
        
        read_localIndices_temp=read_localIndices[i];
        
        coordinate=coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ );
        gdist=coordTransform.edgeDistance(coordinate, NX, NY, NZ );
	
	if (gdist.min() < BoundaryWidth)
	{
		if (gdist.x < BoundaryWidth) {
			/* left boundary */
			if (coordinate.x < BoundaryWidth){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.x)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else a_temp=0.0;
			write_a_x[i]=a_temp;
			write_b_x[i]= b_temp;
			write_k_x[i]=1/k_temp;


			/* right boundary 
			 *starts with half point -> first point shiftet +1 to th right*/
			}else if (gdist.x < BoundaryWidth-1){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.x-1)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else a_temp=0.0;
			write_a_x[i]=a_temp;
			write_b_x[i]= b_temp;
			write_k_x[i]=1/k_temp;
			
			}
			
			
			
			/* half points */
			PositionNorm=(ValueType)(BoundaryWidth-gdist.x-0.5)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ a_temp = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else a_temp=0.0;
			write_a_x_half[i]= a_temp;
			write_b_x_half[i]= b_temp;
			write_k_x_half[i]=1/k_temp;
			
		}
		if (gdist.y < BoundaryWidth) {
			/* left boundary */
			if (coordinate.y < BoundaryWidth){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.y)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_y[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_y[i]=0.0;
			write_b_y[i]= b_temp;
			write_k_y[i]=1/k_temp;
			
			/* right boundary 
			 *starts with half point -> first point shiftet +1 to th right*/
			}else if (gdist.y < BoundaryWidth-1){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.y-1)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_y[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_y[i]=0.0;	
			write_b_y[i]= b_temp;
			write_k_y[i]=1/k_temp;
			
			}
			
			/* half points */
			PositionNorm=(ValueType)(BoundaryWidth-gdist.y-0.5)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_y_half[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_y_half[i]=0.0;
			write_b_y_half[i]= b_temp;
			write_k_y_half[i]=1/k_temp;
		}
		if (gdist.z < BoundaryWidth) {
			/* left boundary */
			if (coordinate.z < BoundaryWidth){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.z)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_z[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_z[i]=0.0;
			write_b_z[i]= b_temp;
			write_k_z[i]=1/k_temp;
			
			/* right boundary 
			 *starts with half point -> first point shiftet +1 to th right*/
			}else if (gdist.z < BoundaryWidth-1){
			PositionNorm=(ValueType)(BoundaryWidth-gdist.z-1)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_z[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_z[i]=0.0;	
			write_b_z[i]= b_temp;
			write_k_z[i]=1/k_temp;
			
			}
			
			/* half points */
			PositionNorm=(ValueType)(BoundaryWidth-gdist.z-0.5)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_z_half[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_z_half[i]=0.0;
			write_b_z_half[i]= b_temp;
			write_k_z_half[i]=1/k_temp;
		}
	
      }
      }
//     

//     /* Release all read and write access */
    read_localIndices.release();
   write_k_x.release();
   write_k_y.release();
   write_k_z.release();
   write_a_x.release();
   write_a_y.release();
   write_a_z.release();
   write_b_x.release();
   write_b_y.release();
   write_b_z.release();
   write_k_x_half.release();
   write_k_y_half.release();
   write_k_z_half.release();
   write_a_x_half.release();
   write_a_y_half.release();
   write_a_z_half.release();
   write_b_x_half.release();
   write_b_y_half.release();
   write_b_z_half.release();
    
	k_x.setContextPtr ( ctx );
	k_y.setContextPtr ( ctx );
	k_z.setContextPtr ( ctx );
	b_x.setContextPtr ( ctx );
	b_y.setContextPtr ( ctx );
	b_y.setContextPtr ( ctx );
	a_z.setContextPtr ( ctx );
	a_z.setContextPtr ( ctx );
	a_z.setContextPtr ( ctx );
	
	k_x_half.setContextPtr ( ctx );
	k_y_half.setContextPtr ( ctx );
	k_z_half.setContextPtr ( ctx );
	b_x_half.setContextPtr ( ctx );
	b_y_half.setContextPtr ( ctx );
	b_y_half.setContextPtr ( ctx );
	a_z_half.setContextPtr ( ctx );
	a_z_half.setContextPtr ( ctx );
	a_z_half.setContextPtr ( ctx );
	
	HOST_PRINT ( dist->getCommunicatorPtr(), "Finished with initialization of the CPML coefficients!\n\n" );

}



