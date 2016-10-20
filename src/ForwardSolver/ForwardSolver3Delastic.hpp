
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
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm);
            
            Acquisition::Seismogram<ValueType> seismogram; //!< Storage of seismogram data
            
        private:
            
            void gatherSeismograms(Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t);
            void applySource(Acquisition::Sources<ValueType>& sources, Wavefields::Wavefields<ValueType>& wavefield,IndexType NT, IndexType t);
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


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
void KITGPI::ForwardSolver::FD3Delastic<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm){
    
    SCAI_REGION( "timestep" )
    
    /* Get references to required modelparameter */
    lama::DenseVector<ValueType>& inverseDensity=model.getInverseDensity();
    lama::DenseVector<ValueType>& M=model.getM();
    lama::DenseVector<ValueType>& Mu=model.getMu();
    
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
    lama::CSRSparseMatrix<ValueType>& A=derivatives.getA();
    lama::CSRSparseMatrix<ValueType>& B=derivatives.getB();
    lama::CSRSparseMatrix<ValueType>& C=derivatives.getC();
    lama::CSRSparseMatrix<ValueType>& D=derivatives.getD();
    lama::CSRSparseMatrix<ValueType>& E=derivatives.getE();
    lama::CSRSparseMatrix<ValueType>& F=derivatives.getF();
    
    /* Init seismograms */
    seismogram.init(receiver, NT, vX.getContextPtr());
    
    common::unique_ptr<lama::Vector> updatePtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
    common::unique_ptr<lama::Vector> vxxPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vyyPtr( vX.newVector() );
    common::unique_ptr<lama::Vector> vzzPtr( vX.newVector() );
    
    lama::Vector& vxx = *vxxPtr;
    lama::Vector& vyy = *vyyPtr;
    lama::Vector& vzz = *vzzPtr;
    
    /* --------------------------------------- */
    /* Start runtime critical part             */
    /* --------------------------------------- */
    
    for ( IndexType t = 0; t < NT; t++ ){
        
        
        if( t % 100 == 0 && t != 0){
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        /* update velocity */
        update = A * Sxx;
        update += E * Sxy;
        update += F * Sxz;
        update.scale(inverseDensity);
        vX += update;
        
        update = D * Sxy;
        update += B * Syy;
        update += F * Syz;
        update.scale(inverseDensity);
        vY += update;
        
        update = D * Sxz;
        update += E * Syz;
        update += C * Szz;
        update.scale(inverseDensity);
        vZ += update;
        
        /* pressure update */
        update = B * vX;
        update += A * vY;
        update.scale(Mu);
        Sxy += update;
        
        update = C * vX;
        update += A * vZ;
        update.scale(Mu);
        Sxz += update;
        
        update = C * vY;
        update += B * vZ;
        update.scale(Mu);
        Syz += update;
        
        vxx = D * vX;
        vyy = E * vY;
        vzz = F * vZ;
        
        update = vxx;
        update += vyy;
        update += vzz;
        update.scale(M);
        
        vxx.scale(Mu);
        vyy.scale(Mu);
        vzz.scale(Mu);
        
        Sxx += update;
        Sxx += 2 * vxx;
        Syy += update;
        Syy += 2 * vyy;
        Szz += update;
        Szz += 2 * vzz;

        /* Apply source and save seismogram */
        applySource(sources,wavefield,NT,t);
        gatherSeismograms(wavefield,NT,t);

    }
    
    /* --------------------------------------- */
    /* Stop runtime critical part             */
    /* --------------------------------------- */
}
