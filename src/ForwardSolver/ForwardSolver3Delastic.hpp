
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
#include "SourceReceiverImpl/FDTD3Delastic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief 3-D elastic forward solver
        template<typename ValueType>
        class FD3Delastic : public ForwardSolver<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD3Delastic(){};
            
            //! Default destructor
            ~FD3Delastic(){};
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType> const& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT) override;
            
            void prepareBoundaryConditions(Configuration::Configuration const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx) override;
            
        private:
            
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Delastic<ValueType> FreeSurface; //!< Free Surface boundary condition class
            using ForwardSolver<ValueType>::useFreeSurface;
            
            BoundaryCondition::ABS3D<ValueType> DampingBoundary; //!< Damping boundary condition class
            using ForwardSolver<ValueType>::useDampingBoundary;
            
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::prepareBoundaryConditions(Configuration::Configuration const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
    /* Prepare Free Surface */
    if(config.get<IndexType>("FreeSurface")){
        useFreeSurface=true;
        FreeSurface.init(dist,derivatives,config.get<IndexType>("NX"),config.get<IndexType>("NY"),config.get<IndexType>("NZ"),config.get<ValueType>("DT"),config.get<ValueType>("DH"));
    }
    
    /* Prepare Damping Boundary */
    if(config.get<IndexType>("DampingBoundary")){
        useDampingBoundary=true;
        DampingBoundary.init(dist,ctx,config.get<IndexType>("NX"),config.get<IndexType>("NY"),config.get<IndexType>("NZ"),config.get<IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"),useFreeSurface);
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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType> const& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType /*DT*/){
    
    SCAI_REGION( "timestep" )
    
    SCAI_ASSERT_ERROR( NT > 0 , " Number of time steps has to be greater than zero. ");
    
    /* Get references to required modelparameter */
    lama::Vector const& inverseDensity=model.getInverseDensity();
    lama::Vector const& pWaveModulus=model.getPWaveModulus();
    lama::Vector const& sWaveModulus=model.getSWaveModulus();
    lama::Vector const& inverseDensityAverageX=model.getInverseDensityAverageX();
    lama::Vector const& inverseDensityAverageY=model.getInverseDensityAverageY();
    lama::Vector const& inverseDensityAverageZ=model.getInverseDensityAverageZ();
    lama::Vector const& sWaveModulusAverageXY=model.getSWaveModulusAverageXY();
    lama::Vector const& sWaveModulusAverageXZ=model.getSWaveModulusAverageXZ();
    lama::Vector const& sWaveModulusAverageYZ=model.getSWaveModulusAverageYZ();
    
    /* Get references to required wavefields */
    lama::Vector & vX=wavefield.getVX();
    lama::Vector & vY=wavefield.getVY();
    lama::Vector & vZ=wavefield.getVZ();
    
    lama::Vector & Sxx=wavefield.getSxx();
    lama::Vector & Syy=wavefield.getSyy();
    lama::Vector & Szz=wavefield.getSzz();
    
    lama::Vector & Syz=wavefield.getSyz();
    lama::Vector & Sxz=wavefield.getSxz();
    lama::Vector & Sxy=wavefield.getSxy();
    
    /* Get references to required derivatives matrixes */
    lama::Matrix const& Dxf=derivatives.getDxf();
    lama::Matrix const& Dzf=derivatives.getDzf();
    lama::Matrix const& Dxb=derivatives.getDxb();
    lama::Matrix const& Dzb=derivatives.getDzb();
    
    lama::Matrix const& DybPressure=derivatives.getDybPressure();
    lama::Matrix const& DybVelocity=derivatives.getDybVelocity();
    lama::Matrix const& DyfPressure=derivatives.getDyfPressure();
    lama::Matrix const& DyfVelocity=derivatives.getDyfVelocity();
    
    SourceReceiverImpl::FDTD3Delastic<ValueType> SourceReceiver(sources,receiver,wavefield);
    
    common::unique_ptr<lama::Vector> updatePtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
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
        update += DybVelocity * Sxy;
        update += Dzb * Sxz;
        vX += update.scale(inverseDensityAverageX);
        
        update = Dxb * Sxy;
        update += DyfVelocity * Syy;
        update += Dzb * Syz;
        vY += update.scale(inverseDensityAverageY);
        
        update = Dxb * Sxz;
        update += DybVelocity * Syz;
        update += Dzf * Szz;
        vZ += update.scale(inverseDensityAverageZ);
        
        /* ----------------*/
        /* pressure update */
        /* ----------------*/
        vxx = Dxb * vX;
        vyy = DybPressure * vY;
        vzz = Dzb * vZ;
        
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
        
        update = DyfPressure * vX;
        update += Dxf * vY;
        Sxy += update.scale(sWaveModulusAverageXY);
        
        update = Dzf * vX;
        update += Dxf * vZ;
        Sxz += update.scale(sWaveModulusAverageXZ);
        
        update = Dzf * vY;
        update += DyfPressure * vZ;
        Syz += update.scale(sWaveModulusAverageYZ);
        
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
