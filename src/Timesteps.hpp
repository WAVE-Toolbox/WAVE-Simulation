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
               lama::Scalar v_factor, lama::Scalar p_factor,
               IndexType NT, lama::Scalar DH_INV, dmemo::CommunicatorPtr comm, dmemo::DistributionPtr /*dist*/ )
{
    SCAI_REGION( "timestep" )
    
    // Invert Density Values before the time stepping
    model.density.invert();
    
    ValueType start_t, end_t;
    ValueType time_sum=0;
    
    // create new Vector(Pointer) with same configuration as vZ
    common::unique_ptr<lama::Vector> updatePtr( wavefield.vX.newVector() );
    // get Reference of VectorPointer
    lama::Vector& update = *updatePtr;
    
    for ( IndexType t = 0; t < NT; t++ )
    {
        if( t % 100 == 0 && t != 0)
        {
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        // update velocity, v_factor is 'DT / DH'
        // velocity z: vZ = vZ + DT / ( DH * rho ) * A * p;
        update=v_factor * A * wavefield.p; // Update=DT / ( DH) * A * p
        wavefield.vZ += update.scale(model.density); // Update+1/RHO
        
        // velocity x: vX = vX + DT / ( DH * rho ) * B * p;
        update=v_factor * B * wavefield.p; // Update=DT / ( DH) * B * p
        wavefield.vX += update.scale(model.density); // Update+1/RHO
        
        // velocity y: vY = vY + DT / ( DH * rho ) * C * p;
        update=v_factor * C * wavefield.p; // Update=DT / ( DH) * C * p
        wavefield.vY += update.scale(model.density); // Update+1/RHO
        
        
        // pressure update
        update =  DH_INV * D * wavefield.vZ;
        update += DH_INV * E * wavefield.vX;
        update += DH_INV * F * wavefield.vY;
        wavefield.p += p_factor * update.scale(model.pi); // p= DT (p_factor) * Update * Model (M)
        
        start_t = common::Walltime::get();
        sources.applySourceLocal(wavefield.p,t,NT);
        receiver.saveSeismogramsLocal(wavefield.p,t,NT);
        end_t = common::Walltime::get();
        time_sum+=(end_t-start_t);
        
    }
    
    HOST_PRINT( comm, "Total time acquisition: " << time_sum << " sec.\n" );
    
}


