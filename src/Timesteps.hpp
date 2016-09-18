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

#include "Modelparameter.hpp"

#include "Wavefields/Wavefields3Dacoustic.hpp"

/*
 *  routine doing NT time steps updating vX, vY, vZ, p
 *  with incoming source
 *  storing seismogram data
 */
template <typename ValueType>
void timesteps( lama::DenseVector<ValueType>& seismogram, lama::DenseVector<ValueType>& source, Modelparameter<ValueType>& model, Wavefields3Dacoustic<ValueType>& wavefield,
               lama::Matrix& A, lama::Matrix& B, lama::Matrix& C, lama::Matrix& D, lama::Matrix& E, lama::Matrix& F,
               lama::Scalar v_factor, lama::Scalar p_factor,
               IndexType NT, lama::Scalar DH_INV, IndexType source_index, IndexType seismogram_index,
               dmemo::CommunicatorPtr comm, dmemo::DistributionPtr /*dist*/ )
{
    SCAI_REGION( "timestep" )
    
    // Invert Density Values before the time stepping
    model.density.invert();
    
    
    // create new Vector(Pointer) with same configuration as vZ
    common::unique_ptr<lama::Vector> helpPtr( wavefield.vX.newVector() );
    // get Reference of VectorPointer
    lama::Vector& help = *helpPtr;
    
    for ( IndexType t = 0; t < NT; t++ )
    {
        if( t % 100 == 0 && t != 0)
        {
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        // update velocity, v_factor is 'DT / DH'
        // velocity z: vZ = vZ + DT / ( DH * rho ) * A * p;
        help=v_factor * A * wavefield.p; // Update=DT / ( DH) * A * p
        wavefield.vZ += help.scale(model.density); // Update+1/RHO
        
        // velocity x: vX = vX + DT / ( DH * rho ) * B * p;
        help=v_factor * B * wavefield.p; // Update=DT / ( DH) * B * p
        wavefield.vX += help.scale(model.density); // Update+1/RHO
        
        // velocity y: vY = vY + DT / ( DH * rho ) * C * p;
        help=v_factor * C * wavefield.p; // Update=DT / ( DH) * C * p
        wavefield.vY += help.scale(model.density); // Update+1/RHO
        
        
        // pressure update
        help =  DH_INV * D * wavefield.vZ;
        help += DH_INV * E * wavefield.vX;
        help += DH_INV * F * wavefield.vY;
        wavefield.p += p_factor * help.scale(model.pi); // p= DT (p_factor) * Update * Model (M)
        
        // update seismogram and pressure with source terms
        // CAUTION: elementwise access by setVal and getVal cause performace issues executed on CUDA
        //          should be used rarely
        // TODO: can do this by index operator[] --> no need for DenseVector<>, can use Vector instead
        wavefield.p.setValue( source_index, wavefield.p.getValue( source_index ) + source.getValue( t ) );
        
        seismogram.setValue( t, wavefield.p.getValue( seismogram_index ) );
        
    }
    
}


