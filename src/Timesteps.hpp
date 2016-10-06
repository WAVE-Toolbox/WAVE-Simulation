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

#include <iostream>

#include "Modelparameter/Modelparameter3Dacoustic.hpp"

#include "Wavefields/Wavefields3Dacoustic.hpp"

#include "Sources.hpp"
#include "Receivers.hpp"

/*
 *  routine doing NT time steps updating vX, vY, vZ, p
 *  with incoming source
 *  storing seismogram data
 */
template <typename ValueType>
void timesteps( Receivers<ValueType>& receiver, Sources<ValueType>& sources, Modelparameter3Dacoustic<ValueType>& model, Wavefields3Dacoustic<ValueType>& wavefield,
               lama::Matrix& A, lama::Matrix& B, lama::Matrix& C, lama::Matrix& D, lama::Matrix& E, lama::Matrix& F,
               IndexType NT, dmemo::CommunicatorPtr comm, dmemo::DistributionPtr /*dist*/ )
{
    
    SCAI_REGION( "timestep" )
    
    model.density.invert(); // Invert Density Values before the time stepping
    
    common::unique_ptr<lama::Vector> updatePtr( wavefield.vX.newVector() ); // create new Vector(Pointer) with same configuration as vZ
    lama::Vector& update = *updatePtr; // get Reference of VectorPointer
    
    for ( IndexType t = 0; t < NT; t++ ){
        
        
        if( t % 100 == 0 && t != 0){
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        
        /* update velocity */
        update=  A * wavefield.p;
        wavefield.vZ += update.scale(model.density);
        
        update= B * wavefield.p;
        wavefield.vX += update.scale(model.density);
        
        update= C * wavefield.p;
        wavefield.vY += update.scale(model.density);
        
        
        /* pressure update */
        update  =  D * wavefield.vZ;
        update +=  E * wavefield.vX;
        update +=  F * wavefield.vY;
        wavefield.p += update.scale(model.pi);
        
        
        /* Apply source and save seismogram */
        sources.applySourceLocal(wavefield.p,t,NT);
        receiver.saveSeismogramsLocal(wavefield.p,t,NT);

        
    }
    
}


