
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"

#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "../Derivatives/Derivatives.hpp"

template<typename ValueType>
class ForwardSolver3Dacoustic : public ForwardSolver<ValueType>
{
public:
    
    ForwardSolver3Dacoustic(){};
    ~ForwardSolver3Dacoustic(){};
    
    void run(Receivers<ValueType>& receiver, Sources<ValueType>& sources, Modelparameter<ValueType>& model, Wavefields<ValueType>& wavefield, Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm);
    
};

template<typename ValueType>
void ForwardSolver3Dacoustic<ValueType>::run(Receivers<ValueType>& receiver, Sources<ValueType>& sources, Modelparameter<ValueType>& model, Wavefields<ValueType>& wavefield, Derivatives<ValueType>& derivatives, IndexType NT, dmemo::CommunicatorPtr comm){
    
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
