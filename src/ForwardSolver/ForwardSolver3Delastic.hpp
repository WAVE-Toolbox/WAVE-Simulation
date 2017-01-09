
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <iostream>

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"


#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "Derivatives/Derivatives.hpp"
#include "BoundaryCondition/FreeSurface3Delastic.hpp"
#include "BoundaryCondition/ABS3D.hpp"
#include "BoundaryCondition/CPML3D.hpp"
#include "SourceReceiverImpl/FDTD3Delastic.hpp"
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
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT) override;
            
            void prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx) override;
            
            
            
        private:
            
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;
            
            BoundaryCondition::ABS3D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;
	    
	    BoundaryCondition::CPML3D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;
            
            
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
    /* Prepare Free Surface */
    if(config.getFreeSurface()){
        useFreeSurface=true;
        FreeSurface.init(dist,derivatives,config.getNX(),config.getNY(),config.getNZ(),config.getDT(),config.getDH());
    }
    
    /* Prepare Damping Boundary */
    if(config.getDampingBoundary()==1){
        useDampingBoundary=true;
        DampingBoundary.init(dist,ctx,config.getNX(),config.getNY(),config.getNZ(),config.getBoundaryWidth(), config.getDampingCoeff(),useFreeSurface);
    }
    if(config.getDampingBoundary()==2){
	useConvPML=true;
	ConvPML.init(dist,ctx,config.getNX(),config.getNY(),config.getNZ(),config.getDT(),config.getDH(),config.getBoundaryWidth(),useFreeSurface,config.getPMLVariables());
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
 \param DT Temporal Sampling intervall in seconds
 */
template<typename ValueType>
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType /*DT*/){   
    SCAI_REGION( "timestep" )
    
    SCAI_ASSERT_ERROR( NT > 0 , " Number of time steps has to be greater than zero. ");
    
    /* Get references to required modelparameter */
    lama::DenseVector<ValueType>const& inverseDensity=model.getInverseDensity();
    lama::DenseVector<ValueType>const& pWaveModulus=model.getPWaveModulus();
    lama::DenseVector<ValueType>const& sWaveModulus=model.getSWaveModulus();
    
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
    lama::CSRSparseMatrix<ValueType>const& Dxf=derivatives.getDxf();
    lama::CSRSparseMatrix<ValueType>const& Dzf=derivatives.getDzf();
    lama::CSRSparseMatrix<ValueType>const& Dxb=derivatives.getDxb();
    lama::CSRSparseMatrix<ValueType>const& Dzb=derivatives.getDzb();
    
    lama::CSRSparseMatrix<ValueType>const& DybPressure=derivatives.getDybPressure();
    lama::CSRSparseMatrix<ValueType>const& DybVelocity=derivatives.getDybVelocity();
    lama::CSRSparseMatrix<ValueType>const& DyfPressure=derivatives.getDyfPressure();
    lama::CSRSparseMatrix<ValueType>const& DyfVelocity=derivatives.getDyfVelocity();
    
    SourceReceiverImpl::FDTD3Delastic<ValueType> SourceReceiver(sources,receiver,wavefield);
    
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
    
    if(useFreeSurface){
        FreeSurface.setModelparameter(model);
    }
    
    dmemo::CommunicatorPtr comm=inverseDensity.getDistributionPtr()->getCommunicatorPtr();

    
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
	if(useConvPML)  {ConvPML.apply_sxx_x(update);}
	
	update_temp = DybVelocity * Sxy;
	if(useConvPML) { ConvPML.apply_sxy_y(update_temp);}
	update += update_temp;
	
        update_temp = Dzb * Sxz;
	if(useConvPML)  { ConvPML.apply_sxz_z(update_temp); }
	update += update_temp;
	
        vX += update.scale(inverseDensity);
        
        update = Dxb * Sxy;
	if(useConvPML) { ConvPML.apply_sxy_x(update);}
	
        update_temp = DyfVelocity * Syy;
	if(useConvPML) { ConvPML.apply_syy_y(update_temp);}
	update += update_temp;
	
        update_temp = Dzb * Syz;
	if(useConvPML) { ConvPML.apply_syz_z(update_temp);}
	update += update_temp;
	
        vY += update.scale(inverseDensity);
        
        update = Dxb * Sxz;
	if(useConvPML) { ConvPML.apply_sxz_x(update);}
	
        update_temp = DybVelocity * Syz;
	if(useConvPML) { ConvPML.apply_syz_y(update_temp);}
	update += update_temp;
	
        update_temp = Dzf * Szz;
	if(useConvPML) { ConvPML.apply_szz_z(update_temp);}
	update += update_temp;
	
        vZ += update.scale(inverseDensity);
        
        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        vzz = Dzb * vZ;
	if(useConvPML) { 
		ConvPML.apply_vxx(vxx);
		ConvPML.apply_vyy(vyy);
		ConvPML.apply_vzz(vzz);
	}
        
        update = vxx;
        update += vyy;
        update += vzz;
        update.scale(pWaveModulus);
        
        Sxx += update;
        Syy += update;
        Szz += update;
        
        update=vyy+vzz;
        Sxx -= 2.0 * update.scale(sWaveModulus);
        update=vxx+vzz;
        Syy -= 2.0 * update.scale(sWaveModulus);
        update=vxx+vyy;
        Szz -= 2.0 * update.scale(sWaveModulus);
        
	//================================
        update = DyfPressure * vX;
	if(useConvPML) { ConvPML.apply_vxy(update);}
	
        update_temp = Dxf * vY;
	if(useConvPML)  {ConvPML.apply_vyx(update_temp);}
	update += update_temp;
	
        Sxy += update.scale(sWaveModulus);
        //====================================
        update = Dzf * vX;
	if(useConvPML) { ConvPML.apply_vxz(update);}
        
        update_temp = Dxf * vZ;
	if(useConvPML)  {ConvPML.apply_vzx(update_temp);}
	update += update_temp;
	
        Sxz += update.scale(sWaveModulus);
        //=========================================
        update = Dzf * vY;
	if(useConvPML) { ConvPML.apply_vyz(update);}
        
	
        update_temp = DyfPressure * vZ;
	if(useConvPML) { ConvPML.apply_vzy(update_temp);}
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
        SourceReceiver.applySource(t);
        SourceReceiver.gatherSeismogram(t);
        
    }
    ValueType end_t = common::Walltime::get();
    HOST_PRINT( comm, "Finished time stepping in " << end_t - start_t << " sec.\n\n" );
    
    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}
