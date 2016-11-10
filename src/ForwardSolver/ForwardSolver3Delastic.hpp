
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <iostream>

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"
#include "../Acquisition/Seismogram.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "Derivatives/Derivatives.hpp"
#include "BoundaryCondition/FreeSurface3Delastic.hpp"
#include "BoundaryCondition/ABS3D.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief 3-D elastic forward solver
        template<typename ValueType>
        class FD3Delastic : public ForwardSolver<ValueType>
        {
            
        public:
            
            /* Default constructor */
            FD3Delastic(){};
            
            /* Default destructor */
            ~FD3Delastic(){};
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm,dmemo::DistributionPtr dist);
            
            void prepareBoundaryConditions(Configuration::Configuration<ValueType> config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx);
            
            using ForwardSolver<ValueType>::seismogram;
            
        private:
            
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;
            
            BoundaryCondition::ABS3D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;
            
            void gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t);
            void applySource(Acquisition::Sources<ValueType>& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t);
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Initialitation of the boundary conditions
 *
 *
 \param config Configuration
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 */
template<typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration<ValueType> config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
    /* Prepare Free Surface */
    if(config.getFreeSurface()){
        useFreeSurface=true;
        FreeSurface.init(dist,derivatives,config.getNX(),config.getNY(),config.getNZ(),config.getDT(),config.getDH());
    }
    
    /* Prepare Damping Boundary */
    if(config.getDampingBoundary()){
        useDampingBoundary=true;
        DampingBoundary.init(dist,ctx,config.getNX(),config.getNY(),config.getNZ(),config.getBoundaryWidth(), config.getDampingCoeff(),useFreeSurface);
    }
    
}

/*! \brief Appling the sources to the wavefield
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param sources Sources to apply
 \param wavefield Wavefields
 \param NT Total number of time steps
 \param t Current time step
 */
template<typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::applySource(Acquisition::Sources<ValueType>& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)
{
    
    IndexType numSourcesLocal=sources.getNumSourcesLocal();
    
    if(numSourcesLocal>0){
        
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& vX=wavefield.getVX();
        lama::DenseVector<ValueType>& vY=wavefield.getVY();
        lama::DenseVector<ValueType>& vZ=wavefield.getVZ();
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
        lama::DenseVector<ValueType>& Szz=wavefield.getSzz();
        
        /* Get reference to sourcesignal storing seismogram */
        Acquisition::Seismogram<ValueType>& signals=sources.getSignals();
        
        /* Get reference to source type of sources */
        lama::DenseVector<ValueType>& SourceType=signals.getTraceType();
        utilskernel::LArray<ValueType>* SourceType_LA=&SourceType.getLocalValues();
        hmemo::WriteAccess<ValueType> read_SourceType_LA(*SourceType_LA);
        
        /* Get reference to coordinates of sources */
        lama::DenseVector<ValueType>& coordinates=signals.getCoordinates();
        utilskernel::LArray<ValueType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<ValueType> read_coordinates_LA(*coordinates_LA);
        
        /* Get reference to storage of source signals */
        lama::DenseMatrix<ValueType>& sourcesSignals=signals.getData();
        lama::DenseStorage<ValueType>* sourcesSignals_DS=&sourcesSignals.getLocalStorage();
        hmemo::HArray<ValueType>* sourcesSignals_HA=&sourcesSignals_DS->getData();
        hmemo::ReadAccess<ValueType> read_sourcesSignals_HA(*sourcesSignals_HA);
        
        /* Get the distribution of the wavefield*/
        dmemo::DistributionPtr dist_wavefield=vX.getDistributionPtr();
        
        IndexType coordinate_global;
        IndexType coordinate_local;
        
        for(IndexType i=0; i<numSourcesLocal; i++){
            coordinate_global=read_coordinates_LA[i];
            coordinate_local=dist_wavefield->global2local(coordinate_global);
            
            switch (IndexType(read_SourceType_LA[i])) {
                case 1:
                    Sxx.getLocalValues()[coordinate_local] = Sxx.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    Syy.getLocalValues()[coordinate_local] = Syy.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    Szz.getLocalValues()[coordinate_local] = Szz.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                case 2:
                    vX.getLocalValues()[coordinate_local] = vX.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                case 3:
                    vY.getLocalValues()[coordinate_local] = vY.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                case 4:
                    vZ.getLocalValues()[coordinate_local] = vZ.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                default:
                    COMMON_THROWEXCEPTION("Source type is unkown")
                    break;
            }
        }
        
        read_coordinates_LA.release();
        read_sourcesSignals_HA.release();
    }
}



/*! \brief Saving seismograms during time stepping
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 \param wavefield Wavefields
 \param NT Total number of time steps
 \param t Current time step
 */
template<typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)
{
    
    IndexType numTracesLocal=seismogram.getNumTracesLocal();
    
    if(numTracesLocal>0){
        
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& vX=wavefield.getVX();
        lama::DenseVector<ValueType>& vY=wavefield.getVY();
        lama::DenseVector<ValueType>& vZ=wavefield.getVZ();
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
        lama::DenseVector<ValueType>& Szz=wavefield.getSzz();
        
        /* Get reference to receiver type of seismogram traces */
        lama::DenseVector<ValueType>& ReceiverType=seismogram.getTraceType();
        utilskernel::LArray<ValueType>* ReceiverType_LA=&ReceiverType.getLocalValues();
        hmemo::WriteAccess<ValueType> read_ReceiverType_LA(*ReceiverType_LA);
        
        /* Get reference to coordinates of seismogram traces */
        lama::DenseVector<ValueType>& coordinates=seismogram.getCoordinates();
        utilskernel::LArray<ValueType>* coordinates_LA=&coordinates.getLocalValues();
        hmemo::WriteAccess<ValueType> read_coordinates_LA(*coordinates_LA);
        
        /* Get reference to storage of seismogram traces */
        lama::DenseMatrix<ValueType>& seismogramData=seismogram.getData();
        lama::DenseStorage<ValueType>* seismogram_DS=&seismogramData.getLocalStorage();
        hmemo::HArray<ValueType>* seismogram_HA=&seismogram_DS->getData();
        hmemo::WriteAccess<ValueType> write_seismogram_HA(*seismogram_HA);
        
        /* Get the distribution of the wavefield*/
        dmemo::DistributionPtr dist_wavefield=vX.getDistributionPtr();
        
        IndexType coordinate_global;
        IndexType coordinate_local;
        ValueType temp;
        
        for(IndexType i=0; i<numTracesLocal; i++){
            coordinate_global=read_coordinates_LA[i];
            coordinate_local=dist_wavefield->global2local(coordinate_global);
            
            switch (IndexType(read_ReceiverType_LA[i])) {
                case 1:
                    temp=Sxx.getLocalValues()[coordinate_local]+Syy.getLocalValues()[coordinate_local]+Szz.getLocalValues()[coordinate_local];
                    write_seismogram_HA[t+NT*i]=temp;
                    break;
                case 2:
                    write_seismogram_HA[t+NT*i]=vX.getLocalValues()[coordinate_local];
                    break;
                case 3:
                    write_seismogram_HA[t+NT*i]=vY.getLocalValues()[coordinate_local];
                    break;
                case 4:
                    write_seismogram_HA[t+NT*i]=vZ.getLocalValues()[coordinate_local];
                    break;
                default:
                    COMMON_THROWEXCEPTION("Receiver type is unkown")
                    break;
            }
        }
        
        read_coordinates_LA.release();
        write_seismogram_HA.release();
    }
}


/*! \brief Running the 3-D elastic foward solver
 *
 * Start the 3-D forward solver as defined by the given parameters
 *
 \param receiver Configuration of the receivers
 \param sources Configuration of the sources
 \param model Configuration of the modelparameter
 \param wavefield Wavefields for the modelling
 \param derivatives Derivations matrices to calculate the spatial derivatives
 \param NT Total number of time steps
 \param comm Communicator
 */
template<typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm,dmemo::DistributionPtr dist){
    
    SCAI_REGION( "timestep" )
    
    /* Get references to required modelparameter */
    lama::DenseVector<ValueType>& inverseDensity=model.getInverseDensity();
    lama::DenseVector<ValueType>& pWaveModulus=model.getPWaveModulus();
    lama::DenseVector<ValueType>& sWaveModulus=model.getSWaveModulus();
    
    /* Get references to required wavefields */
    lama::DenseVector<ValueType>& vX=wavefield.getVX();
    lama::DenseVector<ValueType>& vY=wavefield.getVY();
    lama::DenseVector<ValueType>& vZ=wavefield.getVZ();
    
    lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
    lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
    lama::DenseVector<ValueType>& Szz=wavefield.getSzz();
    
    lama::DenseVector<ValueType>& Syz=wavefield.getSyz();
    lama::DenseVector<ValueType>& Sxz=wavefield.getSxz();
    lama::DenseVector<ValueType>& Sxy=wavefield.getSxy();
    
    /* Get references to required derivatives matrixes */
    lama::CSRSparseMatrix<ValueType>& Dxf=derivatives.getDxf();
    lama::CSRSparseMatrix<ValueType>& Dzf=derivatives.getDzf();
    lama::CSRSparseMatrix<ValueType>& Dxb=derivatives.getDxb();
    lama::CSRSparseMatrix<ValueType>& Dzb=derivatives.getDzb();
    
    lama::CSRSparseMatrix<ValueType>& DybPressure=derivatives.getDybPressure();
    lama::CSRSparseMatrix<ValueType>& DybVelocity=derivatives.getDybVelocity();
    lama::CSRSparseMatrix<ValueType>& DyfPressure=derivatives.getDyfPressure();
    lama::CSRSparseMatrix<ValueType>& DyfVelocity=derivatives.getDyfVelocity();
    
    /* Init seismograms */
    seismogram.init(receiver, NT, vX.getContextPtr());
    
    common::unique_ptr<lama::Vector> updatePtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
    common::unique_ptr<lama::Vector> update_tempPtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update_temp = *update_tempPtr; // get Reference of VectorPointer
    
    common::unique_ptr<lama::Vector> vxxPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vzzPtr( vX.newVector() );
    
    lama::Vector& vxx = *vxxPtr;
    lama::Vector& vyy = *vyyPtr;
    lama::Vector& vzz = *vzzPtr;
    

    
    
    lama::DenseVector<ValueType> psi_vxx(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vyx(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vzx(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vxy(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vyy(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vzy(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vxz(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vyz(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_vzz(vX.getDistributionPtr());

    
    lama::DenseVector<ValueType> psi_sxx_x(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_sxy_x(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_sxz_x(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_sxy_y(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_syy_y(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_syz_y(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_sxz_z(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_syz_z(vX.getDistributionPtr());
    lama::DenseVector<ValueType> psi_szz_z(vX.getDistributionPtr());

    
    
    if(useFreeSurface){
        FreeSurface.setModelparameter(model);
    }
    
    
    /* --------------------------------------- */
    /*PML init			           */
    /* --------------------------------------- */
    IndexType BoundaryWidth=10;
    IndexType DH=50;
    
    
    IndexType NPower=4;
    ValueType PI=3.14159265359;
    ValueType K_Max_Pml=1;
    ValueType FPml=5.0;
    ValueType RCoef=0.0008;
    ValueType VPml=3500;
    
    ValueType alpha_max_Pml=2.0 * PI * (FPml/2.0);
    
    ValueType d0;
    ValueType PositionNorm;
    
   d0 = - (NPower + 1) * VPml * log(RCoef) / (2.0 * BoundaryWidth*DH);
    

  lama::DenseVector<ValueType> k_x(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_x_LA=&k_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_x(*k_x_LA);
   
  lama::DenseVector<ValueType> b_x(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_x_LA=&b_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_x(*b_x_LA);
   
  lama::DenseVector<ValueType> a_x(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_x_LA=&a_x.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_x(*a_x_LA);
  
  ValueType k_temp=0.0;
  ValueType b_temp=0.0;
  ValueType a_temp=0.0;
  
  lama::DenseVector<ValueType> k_y(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_y_LA=&k_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_y(*k_y_LA);
   
  lama::DenseVector<ValueType> b_y(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_y_LA=&b_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_y(*b_y_LA);
   
  lama::DenseVector<ValueType> a_y(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_y_LA=&a_y.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_y(*a_y_LA);
  
   
  lama::DenseVector<ValueType> k_z(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_z_LA=&k_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_z(*k_z_LA);
   
  lama::DenseVector<ValueType> b_z(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_z_LA=&b_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_z(*b_z_LA);
   
  lama::DenseVector<ValueType> a_z(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_z_LA=&a_z.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_z(*a_z_LA);
  
  
  
    lama::DenseVector<ValueType> k_x_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_x_half_LA=&k_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_x_half(*k_x_half_LA);
   
  lama::DenseVector<ValueType> b_x_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_x_half_LA=&b_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_x_half(*b_x_half_LA);
   
  lama::DenseVector<ValueType> a_x_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_x_half_LA=&a_x_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_x_half(*a_x_half_LA);
   
  lama::DenseVector<ValueType> k_y_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_y_half_LA=&k_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_y_half(*k_y_half_LA);
   
  lama::DenseVector<ValueType> b_y_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_y_half_LA=&b_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_y_half(*b_y_half_LA);
   
  lama::DenseVector<ValueType> a_y_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_y_half_LA=&a_y_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_y_half(*a_y_half_LA);
   
  lama::DenseVector<ValueType> k_z_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* k_z_half_LA=&k_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_k_z_half(*k_z_half_LA);
   
  lama::DenseVector<ValueType> b_z_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* b_z_half_LA=&b_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_b_z_half(*b_z_half_LA);
   
  lama::DenseVector<ValueType> a_z_half(vX.getDistributionPtr());
  utilskernel::LArray<ValueType>* a_z_half_LA=&a_z_half.getLocalValues();
  hmemo::WriteAccess<ValueType> write_a_z_half(*a_z_half_LA);
  
 
   
  ValueType alpha_prime=0.0;
  ValueType d=0.0;
  
   
   
   
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices);
    
    IndexType numLocalIndices=localIndices.size(); // Number of local indices
    
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp; // Temporary storage, so we do not have to access the array
   
    Acquisition::Coordinates<ValueType> coordTransform;
    Acquisition::coordinate3D coordinate;
    Acquisition::coordinate3D gdist;
    
    //IndexType distgrid;
    IndexType NX=100, NY=100, NZ=100;
    ValueType DT=2.0e-03;
    

    
       for( IndexType i=0; i<numLocalIndices; i++ ) {
        
        read_localIndices_temp=read_localIndices[i];
        
        coordinate=coordTransform.index2coordinate(read_localIndices_temp, NX, NY, NZ );
        gdist=coordTransform.edgeDistance(coordinate, NX, NY, NZ );
        
	write_k_x[i]=1.0;
	write_b_x[i]=1.0;
	write_a_x[i]=1.0;
	write_k_y[i]=1.0;
	write_b_y[i]=1.0;
	write_a_y[i]=1.0;
	write_k_z[i]=1.0;
	write_b_z[i]=1.0;
	write_a_z[i]=1.0;
	write_k_x_half[i]=1.0;
	write_b_x_half[i]=1.0;
	write_a_x_half[i]=1.0;
	write_k_y_half[i]=1.0;
	write_b_y_half[i]=1.0;
	write_a_y_half[i]=1.0;
	write_k_z_half[i]=1.0;
	write_b_z_half[i]=1.0;
	write_a_z_half[i]=1.0;
	
	
	
	if (gdist.min() < BoundaryWidth)
	{
		if (gdist.x < BoundaryWidth) {
			/* left boundary */
			if (coordinate.x < BoundaryWidth){
			PositionNorm=(BoundaryWidth-gdist.x)/BoundaryWidth;
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

			  if((coordinate.x==9)&&(coordinate.z==0))
			std::cout  << "k:" << k_temp << "  b" << b_temp << "  a" << a_temp << "  d" << d << " alpha " << alpha_prime << "pos " << PositionNorm <<coordinate.x<< gdist.x << std::endl; 
			
			/* right boundary 
			 *starts with half point -> first point shiftet +1 to th right*/
			}else if (gdist.x < BoundaryWidth-1){
			PositionNorm=(BoundaryWidth-gdist.x-1)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_x[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_x[i]=0.0;	
			write_b_x[i]= b_temp;
			write_k_x[i]=1/k_temp;
			
			}
			
			/* half points */
			PositionNorm=(BoundaryWidth-gdist.x-0.5)/BoundaryWidth;
			d = d0 * pow(PositionNorm,NPower);
			k_temp = 1.0 + (K_Max_Pml - 1.0) * pow(PositionNorm,NPower);
			alpha_prime = alpha_max_Pml * (1.0 - PositionNorm);
			b_temp = exp(- (d / k_temp + alpha_prime) * DT);
			/* avoid division by zero outside the PML */
			if(abs(d) > 1.0e-6){ write_a_x_half[i] = d * (b_temp - 1.0) / (k_temp * (d + k_temp * alpha_prime));}
			else write_a_x_half[i]=0.0;
			write_b_x_half[i]= b_temp;
			write_k_x_half[i]=1/k_temp;
		}
		if (gdist.y < BoundaryWidth) {
			/* left boundary */
			if (coordinate.y < BoundaryWidth){
			PositionNorm=(BoundaryWidth-gdist.y)/BoundaryWidth;
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
			PositionNorm=(BoundaryWidth-gdist.y-1)/BoundaryWidth;
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
			PositionNorm=(BoundaryWidth-gdist.y-0.5)/BoundaryWidth;
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
			PositionNorm=(BoundaryWidth-gdist.z)/BoundaryWidth;
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
			PositionNorm=(BoundaryWidth-gdist.z-1)/BoundaryWidth;
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
			PositionNorm=(BoundaryWidth-gdist.z-0.5)/BoundaryWidth;
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
    
    k_z.setContextPtr ( vX.getContextPtr() );
    
    b_x_half.setContextPtr ( vX.getContextPtr() );
    a_x_half.setContextPtr ( vX.getContextPtr() );
    
    /* --------------------------------------- */
    /*PML init end			           */
    /* --------------------------------------- */
  //  std::cout << "my vector x looks like: " << b_x << std::endl;
    
    /* --------------------------------------- */
    /* Start runtime critical part             */
    /* --------------------------------------- */
    
    HOST_PRINT( comm, "Start time stepping\n" );
    ValueType start_t = common::Walltime::get();
    for ( IndexType t = 0; t < NT; t++ ){
        
        
        if( t % 100 == 0 && t != 0){
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        

 
        
        /* ----------------*/
        /* update velocity */
        /* ----------------*/
        update = Dxf * Sxx;
	psi_sxx_x=psi_sxx_x.scale(b_x_half)+update.scale(a_x_half);
	update = update.scale(k_x_half) + psi_sxx_x;
	
	update_temp = DybVelocity * Sxy;
	psi_sxy_y=psi_sxy_y.scale(b_y)+update_temp.scale(a_y);
	update_temp = update.scale(k_y) + psi_sxy_y;
	update += update_temp;
	
        update_temp = Dzb * Sxz;
	psi_sxz_z=psi_sxz_z.scale(b_z)+update_temp.scale(a_z);
	update_temp = update.scale(k_z) + psi_sxz_z;
	update += update_temp;
	
        vX += update.scale(inverseDensity);
        
        
	update = Dxb * Sxy;
	psi_sxy_x=psi_sxy_x.scale(b_x)+update.scale(a_x);
	update = update.scale(k_x) + psi_sxy_x;
	
        update_temp = DyfVelocity * Syy;
	psi_syy_y=psi_syy_y.scale(b_y_half)+update_temp.scale(a_y_half);
	update_temp = update.scale(k_y_half) + psi_syy_y;
	update += update_temp;
	
        update_temp = Dzb * Syz;
	psi_syz_z=psi_syz_z.scale(b_z)+update_temp.scale(a_z);
	update_temp = update.scale(k_z) + psi_syz_z;
	update += update_temp;
	
	vY += update.scale(inverseDensity);
        
	
        update = Dxb * Sxz;
	psi_sxz_x=psi_sxz_x.scale(b_x)+update.scale(a_x);
	update = update.scale(k_x) + psi_sxz_x;
	
        update_temp = DybVelocity * Syz;
	psi_syz_y=psi_syz_y.scale(b_y)+update_temp.scale(a_y);
	update_temp = update.scale(k_y) + psi_syz_y;
	update += update_temp;
	
        update_temp = Dzf * Szz;
	psi_szz_z=psi_szz_z.scale(b_z_half)+update_temp.scale(a_z_half);
	update_temp = update.scale(k_z_half) + psi_szz_z;
	update += update_temp;
	
        vZ += update.scale(inverseDensity);
        
	
        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
	psi_vxx=psi_vxx.scale(b_x) + vxx.scale(a_x);
	vxx=vxx.scale(k_x) + psi_vxx;
        vyy = DybPressure * vY;
	psi_vyy=psi_vyy.scale(b_y) + vyy.scale(a_y);
	vyy=vyy.scale(k_y) + psi_vyy;
        vzz = Dzb * vZ;
	psi_vzz=psi_vzz.scale(b_z) + vzz.scale(a_z);
	vzz=vzz.scale(k_z) + psi_vzz;
        
        update = vxx;
	update += vyy;
	update += vzz;
        update.scale(pWaveModulus);
        
        Sxx += update;
        Syy += update;
        Szz += update;
        
        update=vyy+vzz;
        Sxx -= 2 * update.scale(sWaveModulus);
        update=vxx+vzz;
        Syy -= 2 * update.scale(sWaveModulus);
        update=vxx+vyy;
        Szz -= 2 * update.scale(sWaveModulus);
        
	//================================
        update = DyfPressure * vX;
	psi_vxy=psi_vxy.scale(b_y_half) + update.scale(a_y_half);
	update=update.scale(k_y_half) + psi_vxy;
	
        update_temp = Dxf * vY;
	psi_vyx=psi_vyx.scale(b_x_half)+update_temp.scale(a_x_half);
	update_temp = update.scale(k_x_half) + psi_vyx;
	update += update_temp;
	
        Sxy += update.scale(sWaveModulus);
        //====================================
        update = Dzf * vX;
	psi_vxz=psi_vxz.scale(b_z_half) + update.scale(a_z_half);
	update=update.scale(k_z_half) + psi_vxz;
	
        update_temp = Dxf * vZ;
	psi_vzx=psi_vzx.scale(b_x_half)+update_temp.scale(a_x_half);
	update_temp = update.scale(k_x_half) + psi_vzx;
	update += update_temp;
	
        Sxz += update.scale(sWaveModulus);
        //=========================================
        update = Dzf * vY;
	psi_vyz=psi_vyz.scale(b_z_half) + update.scale(a_z_half);
	update=update.scale(k_z_half) + psi_vyz;
	
        update_temp = DyfPressure * vZ;
	psi_vzy=psi_vzy.scale(b_y_half)+update_temp.scale(a_y_half);
	update_temp = update.scale(k_y_half) + psi_vzy;
	update += update_temp;
	
        Syz += update.scale(sWaveModulus);
        
        /* Apply free surface to stress update */
        if(useFreeSurface){
            update=vxx+vzz;
            FreeSurface.apply(update,Sxx,Syy,Szz);
        }
        
        /* Apply the damping boundary */
        if(useDampingBoundary){
            DampingBoundary.apply(Sxx,Syy,Szz,Sxy,Sxz,Syz,vX,vY,vZ);
        }
        
        /* Apply source and save seismogram */
        applySource(sources,wavefield,NT,t);
        gatherSeismograms(wavefield,NT,t);
        
    }
    ValueType end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );
    
    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}
