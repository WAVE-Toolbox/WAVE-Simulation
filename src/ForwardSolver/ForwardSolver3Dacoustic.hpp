
#pragma once

#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <iostream>

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "../Derivatives/Derivatives.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief 3-D Acoustic forward solver
        template<typename ValueType>
        class FD3Dacoustic : public ForwardSolver<ValueType>
        {
        public:
            
            /* Default constructor */
            FD3Dacoustic(){};
            
            /* Default destructor */
            ~FD3Dacoustic(){};
            
            void run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm);
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */


/*! \brief Running the 3-D acoustic foward solver
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
void KITGPI::ForwardSolver::FD3Dacoustic<ValueType>::run(Acquisition::Receivers<ValueType>& receiver, Acquisition::Sources<ValueType>& sources, Modelparameter::Modelparameter<ValueType>& model, Wavefields::Wavefields<ValueType>& wavefield, Derivatives::Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm){
    
    SCAI_REGION( "timestep" )
    
    /* Get references to required modelparameter */
    lama::DenseVector<ValueType>& inverseDensity=*model.getInverseDensity();
    lama::DenseVector<ValueType>& M=*model.getM();
    
    /* Get references to required wavefields */
    lama::DenseVector<ValueType>& vX=*wavefield.getVX();
    lama::DenseVector<ValueType>& vY=*wavefield.getVY();
    lama::DenseVector<ValueType>& vZ=*wavefield.getVZ();
    lama::DenseVector<ValueType>& p=*wavefield.getP();
    
    /* Get references to required derivatives matrixes */
    lama::CSRSparseMatrix<ValueType>& A=*derivatives.getA();
    lama::CSRSparseMatrix<ValueType>& B=*derivatives.getB();
    lama::CSRSparseMatrix<ValueType>& C=*derivatives.getC();
    lama::CSRSparseMatrix<ValueType>& D=*derivatives.getD();
    lama::CSRSparseMatrix<ValueType>& E=*derivatives.getE();
    lama::CSRSparseMatrix<ValueType>& F=*derivatives.getF();
    
    common::unique_ptr<lama::Vector> updatePtr( vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
    for ( IndexType t = 0; t < NT; t++ ){
        
        
        if( t % 100 == 0 && t != 0){
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        
        /* update velocity */
        update=  A * p;
        vZ += update.scale(inverseDensity);
        
        update= B * p;
        vX += update.scale(inverseDensity);
        
        update= C * p;
        vY += update.scale(inverseDensity);
        
        
        /* pressure update */
        update  =  D * vZ;
        update +=  E * vX;
        update +=  F * vY;
        p += update.scale(M);
        
        
        /* Apply source and save seismogram */
        sources.applySourceLocal(p,t,NT);
        receiver.saveSeismogramsLocal(p,t,NT);
        
    }
}
