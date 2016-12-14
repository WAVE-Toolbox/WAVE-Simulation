
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
#include "BoundaryCondition/FreeSurface3Dvisco.hpp"
#include "BoundaryCondition/ABS3D.hpp"
#include "SourceReceiverImpl/FDTD3Delastic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief 3-D visco forward solver
        template<typename ValueType>
        class FD3Dvisco : public ForwardSolver<ValueType>
        {
            
        public:
            
            /* Default constructor */
            FD3Dvisco(){};
            
            /* Default destructor */
            ~FD3Dvisco(){};
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT) override;
            
            void prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx) override;
            
        private:
            
            /* Boundary Conditions */
            BoundaryCondition::FreeSurface3Dvisco<ValueType> FreeSurface; //!< Free Surface boundary condition class
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
void KITGPI::ForwardSolver::FD3Dvisco<ValueType>::prepareBoundaryConditions(Configuration::Configuration<ValueType> const& config, Derivatives::Derivatives<ValueType>& derivatives,dmemo::DistributionPtr dist, hmemo::ContextPtr ctx){
    
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

/*! \brief Running the 3-D visco-elastic foward solver
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
void KITGPI::ForwardSolver::FD3Dvisco<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType> const& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>const& derivatives, IndexType NT, ValueType DT){
    
    SCAI_REGION( "timestep" )
    
    SCAI_ASSERT_ERROR( NT > 0 , " Number of time steps has to be greater than zero. ");
    
    /* Get references to required modelparameter */
    model.prepareForModelling();
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
    
    lama::DenseVector<ValueType>& Rxx=wavefield.getRxx();
    lama::DenseVector<ValueType>& Ryy=wavefield.getRyy();
    lama::DenseVector<ValueType>& Rzz=wavefield.getRzz();
    lama::DenseVector<ValueType>& Ryz=wavefield.getRyz();
    lama::DenseVector<ValueType>& Rxz=wavefield.getRxz();
    lama::DenseVector<ValueType>& Rxy=wavefield.getRxy();
    
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
   
    lama::DenseVector<ValueType> update2(vX.getDistributionPtr());

    common::unique_ptr<lama::Vector> vxxPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vzzPtr( vX.newVector() );
    
    lama::Vector& vxx = *vxxPtr;
    lama::Vector& vyy = *vyyPtr;
    lama::Vector& vzz = *vzzPtr;
        
    lama::DenseVector<ValueType>const& tauS=model.getTauS();
    lama::DenseVector<ValueType>const& tauP=model.getTauP();
    
    IndexType numRelaxationMechanisms=model.getNumRelaxationMechanisms(); // = Number of relaxation mechanisms
    ValueType relaxationTime=1.0/(2.0*M_PI*model.getRelaxationFrequency()); // = 1 / ( 2 * Pi * f_relax )
    ValueType inverseRelaxationTime=1.0/relaxationTime; // = 1 / relaxationTime
    ValueType viscoCoeff1=(1.0-DT/(2.0*relaxationTime)); // = 1 - DT / ( 2 * tau_Sigma_l )
    ValueType viscoCoeff2=1.0/(1.0+DT/(2.0*relaxationTime)); // = ( 1.0 + DT / ( 2 * tau_Sigma_l ) ) ^ - 1
    ValueType DThalf=DT/2.0; // = DT / 2.0
    
    lama::DenseVector<ValueType> onePlusLtauP(vX.getDistributionPtr()); // = ( 1 + L * tauP )
    lama::DenseVector<ValueType> onePlusLtauS(vX.getDistributionPtr()); // = ( 1 + L * tauS )
    
    onePlusLtauP = 1.0;
    onePlusLtauP += numRelaxationMechanisms * tauP;

    onePlusLtauS = 1.0;
    onePlusLtauS += numRelaxationMechanisms * tauS;
    
    if(useFreeSurface){
        FreeSurface.setModelparameter(model,onePlusLtauP,onePlusLtauS);
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
        vX += update.scale(inverseDensity);
        
        update = Dxb * Sxy;
        update += DyfVelocity * Syy;
        update += Dzb * Syz;
        vY += update.scale(inverseDensity);
        
        update = Dxb * Sxz;
        update += DybVelocity * Syz;
        update += Dzf * Szz;
        vZ += update.scale(inverseDensity);
        
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
        
        update2=inverseRelaxationTime*update;
        update2.scale(tauP);
        
        Sxx += DThalf * Rxx;
        Rxx *= viscoCoeff1;
        Rxx -= update2;
        
        Syy += DThalf * Ryy;
        Ryy *= viscoCoeff1;
        Ryy -= update2;
        
        Szz += DThalf * Rzz;
        Rzz *= viscoCoeff1;
        Rzz -= update2;
        
        update.scale(onePlusLtauP);
        Sxx += update;
        Syy += update;
        Szz += update;
        
        
        /* Update Sxx and Rxx */
        update=vyy+vzz;
        update.scale(sWaveModulus);
        update *= 2.0;

        update2=inverseRelaxationTime* update;
        Rxx += update2.scale(tauS);
        Sxx -= update.scale(onePlusLtauS);

        Rxx *= viscoCoeff2;
        Sxx += DThalf * Rxx;
        
        
        /* Update Syy and Ryy */
        update=vxx+vzz;
        update.scale(sWaveModulus);
        update *= 2.0;

        update2=inverseRelaxationTime* update;
        Ryy += update2.scale(tauS);
        Syy -= update.scale(onePlusLtauS);
        
        Ryy *= viscoCoeff2;
        Syy += DThalf * Ryy;
        
        
        /* Update Szz and Szz */
        update=vxx+vyy;
        update.scale(sWaveModulus);
        update *= 2.0;

        update2=inverseRelaxationTime* update;
        Rzz += update2.scale(tauS);
        Szz -= update.scale(onePlusLtauS);
        
        Rzz *= viscoCoeff2;
        Szz += DThalf * Rzz;
        
        
        /* Update Sxy and Rxy*/
        Sxy += DThalf * Rxy;
        Rxy *= viscoCoeff1;

        update = DyfPressure * vX;
        update += Dxf * vY;
        update.scale(sWaveModulus);

        update2 = inverseRelaxationTime * update;
        Rxy -= update2.scale(tauS);
        Sxy += update.scale(onePlusLtauS);
        
        Rxy *= viscoCoeff2;
        Sxy += DThalf * Rxy;

        
        /* Update Sxz and Rxz */
        Sxz += DThalf * Rxz;
        Rxz *= viscoCoeff1;

        update = Dzf * vX;
        update += Dxf * vZ;
        update.scale(sWaveModulus);
        
        update2 = inverseRelaxationTime * update;
        Rxz -= update2.scale(tauS);
        Sxz += update.scale(onePlusLtauS);

        Rxz *= viscoCoeff2;
        Sxz += DThalf * Rxz;
        
        
        /* Update Syz and Syz */
        Syz += DThalf * Ryz;
        Ryz *= viscoCoeff1;

        update = Dzf * vY;
        update += DyfPressure * vZ;
        update.scale(sWaveModulus);
        
        update2 = inverseRelaxationTime * update;
        Ryz -= update2.scale(tauS);
        Syz += update.scale(onePlusLtauS);

        Ryz *= viscoCoeff2;
        Syz += DThalf * Ryz;
        
        /* Apply free surface to stress update */
        if(useFreeSurface){
            update=vxx+vzz;
            FreeSurface.apply(update,update2,Sxx,Syy,Szz,Rxx,Ryy,Rzz);
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
