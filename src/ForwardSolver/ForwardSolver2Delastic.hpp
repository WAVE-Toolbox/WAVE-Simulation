
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
#include "BoundaryCondition/FreeSurface2Delastic.hpp"
#include "BoundaryCondition/ABS2D.hpp"
#include "BoundaryCondition/CPML2D.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief 3-D elastic forward solver
        template<typename ValueType>
        class FD2Delastic : public ForwardSolver<ValueType>
        {
            
        public:
            
            /* Default constructor */
            FD2Delastic(){};
            
            /* Default destructor */
            ~FD2Delastic(){};
            
            void run(Acquisition::Receivers<ValueType> const& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT) override;
            
            void prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx) override;
            
            using ForwardSolver<ValueType>::seismogram;
            
        private:
            
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface2Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;
            
            BoundaryCondition::ABS2D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;
            
	    BoundaryCondition::CPML2D<ValueType> ConvPML; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useConvPML;

            void gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t) override;
            void applySource(Acquisition::Sources<ValueType> const& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t) override;
            
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::applySource(Acquisition::Sources<ValueType> const& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)
{
    
    IndexType numSourcesLocal=sources.getNumSourcesLocal();
    
    if(numSourcesLocal>0){
        
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& vX=wavefield.getVX();
        lama::DenseVector<ValueType>& vY=wavefield.getVY();
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
        
        /* Get reference to sourcesignal storing seismogram */
        const Acquisition::Seismogram<ValueType>& signals=sources.getSignals();
        
        /* Get reference to source type of sources */
        const lama::DenseVector<IndexType>& SourceType=signals.getTraceType();
        const utilskernel::LArray<IndexType>* SourceType_LA=&SourceType.getLocalValues();
        const hmemo::ReadAccess<IndexType> read_SourceType_LA(*SourceType_LA);
        
        /* Get reference to coordinates of sources */
        const lama::DenseVector<IndexType>& coordinates=signals.getCoordinates();
        const utilskernel::LArray<IndexType>* coordinates_LA=&coordinates.getLocalValues();
        const hmemo::ReadAccess<IndexType> read_coordinates_LA(*coordinates_LA);
        
        /* Get reference to storage of source signals */
        const lama::DenseMatrix<ValueType>& sourcesSignals=signals.getData();
        const lama::DenseStorage<ValueType>* sourcesSignals_DS=&sourcesSignals.getLocalStorage();
        const hmemo::HArray<ValueType>* sourcesSignals_HA=&sourcesSignals_DS->getData();
        const hmemo::ReadAccess<ValueType> read_sourcesSignals_HA(*sourcesSignals_HA);
        
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
                    break;
                case 2:
                    vX.getLocalValues()[coordinate_local] = vX.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                case 3:
                    vY.getLocalValues()[coordinate_local] = vY.getLocalValues()[coordinate_local] + read_sourcesSignals_HA[t+NT*i];
                    break;
                default:
                    COMMON_THROWEXCEPTION("Source type is unkown")
                    break;
            }
        }
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t)
{
    
    IndexType numTracesLocal=seismogram.getNumTracesLocal();
    
    if(numTracesLocal>0){
        
        /* Get reference to wavefields */
        lama::DenseVector<ValueType>& vX=wavefield.getVX();
        lama::DenseVector<ValueType>& vY=wavefield.getVY();
        lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
        lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
        
        /* Get reference to receiver type of seismogram traces */
        const lama::DenseVector<IndexType>& ReceiverType=seismogram.getTraceType();
        const utilskernel::LArray<IndexType>* ReceiverType_LA=&ReceiverType.getLocalValues();
        const hmemo::ReadAccess<IndexType> read_ReceiverType_LA(*ReceiverType_LA);
        
        /* Get reference to coordinates of seismogram traces */
        const lama::DenseVector<IndexType>& coordinates=seismogram.getCoordinates();
        const utilskernel::LArray<IndexType>* coordinates_LA=&coordinates.getLocalValues();
        const hmemo::ReadAccess<IndexType> read_coordinates_LA(*coordinates_LA);
        
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
                    temp=Sxx.getLocalValues()[coordinate_local]+Syy.getLocalValues()[coordinate_local];
                    write_seismogram_HA[t+NT*i]=temp;
                    break;
                case 2:
                    write_seismogram_HA[t+NT*i]=vX.getLocalValues()[coordinate_local];
                    break;
                case 3:
                    write_seismogram_HA[t+NT*i]=vY.getLocalValues()[coordinate_local];
                    break;
                default:
                    COMMON_THROWEXCEPTION("Receiver type is unkown")
                    break;
            }
        }
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
void KITGPI::ForwardSolver::FD2Delastic<ValueType>::run(Acquisition::Receivers<ValueType> const& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType /*DT*/){
    
    SCAI_REGION( "timestep" )
    
    SCAI_ASSERT_ERROR( NT > 0 , " Number of time steps has to be greater than zero. ");
    
    /* Get references to required modelparameter */
    lama::DenseVector<ValueType>const& inverseDensity=model.getInverseDensity();
    lama::DenseVector<ValueType>const& pWaveModulus=model.getPWaveModulus();
    lama::DenseVector<ValueType>const& sWaveModulus=model.getSWaveModulus();
    
    /* Get references to required wavefields */
    lama::DenseVector<ValueType>& vX=wavefield.getVX();
    lama::DenseVector<ValueType>& vY=wavefield.getVY();
    
    lama::DenseVector<ValueType>& Sxx=wavefield.getSxx();
    lama::DenseVector<ValueType>& Syy=wavefield.getSyy();
    
    lama::DenseVector<ValueType>& Sxy=wavefield.getSxy();
    
    /* Get references to required derivatives matrixes */
    lama::CSRSparseMatrix<ValueType>const& Dxf=derivatives.getDxf();
    lama::CSRSparseMatrix<ValueType>const& Dxb=derivatives.getDxb();
    
    lama::CSRSparseMatrix<ValueType>const& DybPressure=derivatives.getDybPressure();
    lama::CSRSparseMatrix<ValueType>const& DybVelocity=derivatives.getDybVelocity();
    lama::CSRSparseMatrix<ValueType>const& DyfPressure=derivatives.getDyfPressure();
    lama::CSRSparseMatrix<ValueType>const& DyfVelocity=derivatives.getDyfVelocity();
    
    /* Init seismograms */
    seismogram.init(receiver, NT, vX.getContextPtr());
    
    common::unique_ptr<lama::Vector> updatePtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
    common::unique_ptr<lama::Vector> update_tempPtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update_temp = *update_tempPtr; // get Reference of VectorPointer
    
    common::unique_ptr<lama::Vector> vxxPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    
    lama::Vector& vxx = *vxxPtr;
    lama::Vector& vyy = *vyyPtr;
    
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
	
        vX += update.scale(inverseDensity);
        
        update = Dxb * Sxy;
	if(useConvPML) { ConvPML.apply_sxy_x(update);}
	
        update_temp = DyfVelocity * Syy;
	if(useConvPML) { ConvPML.apply_syy_y(update_temp);}
	update += update_temp;
	
        vY += update.scale(inverseDensity);
        
        
        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        if(useConvPML) {
	ConvPML.apply_vxx(vxx);
	ConvPML.apply_vyy(vyy);
	}
        
        update = vxx;
        update += vyy;
        update.scale(pWaveModulus);
        
        Sxx += update;
        Syy += update;
        
        Sxx -= 2.0 * vyy.scale(sWaveModulus);
        Syy -= 2.0 * vxx.scale(sWaveModulus);
        
        update = DyfPressure * vX;
	if(useConvPML) { ConvPML.apply_vxy(update);}
	
        update_temp = Dxf * vY;
	if(useConvPML)  {ConvPML.apply_vyx(update_temp);}
	update += update_temp;
	
        Sxy += update.scale(sWaveModulus);
        
        
        /* Apply free surface to stress update */
        if(useFreeSurface){
            FreeSurface.apply(vxx,Sxx,Syy);
        }
        
        /* Apply the damping boundary */
        if(useDampingBoundary){
            DampingBoundary.apply(Sxx,Syy,Sxy,vX,vY);
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
