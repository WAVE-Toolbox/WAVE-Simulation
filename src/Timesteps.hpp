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


/*
 *  routine doing NT time steps updating vX, vY, vZ, p
 *  with incoming source
 *  storing seismogram data
 */
template <typename ValueType>
void timesteps( lama::DenseVector<ValueType>& seismogram, lama::DenseVector<ValueType>& source, lama::DenseVector<ValueType>& p,
               lama::Vector& vX, lama::Vector& vY, lama::Vector& vZ,
               lama::Matrix& A, lama::Matrix& B, lama::Matrix& C, lama::Matrix& D, lama::Matrix& E, lama::Matrix& F,
               lama::Scalar v_factor, lama::Scalar p_factor,
               IndexType NT, lama::Scalar DH_INV, IndexType source_index, IndexType seismogram_index,
               dmemo::CommunicatorPtr comm, dmemo::DistributionPtr /*dist*/ )
{
    SCAI_REGION( "timestep" )
    
    for ( IndexType t = 0; t < NT; t++ )
    {
        if( t % 100 == 0 && t != 0)
        {
            HOST_PRINT( comm, "Calculating time step " << t << " from " << NT << "\n" );
        }
        
        // update velocity, v_factor is 'DT / DH / rho'
        // velocity z: vZ = vZ + DT / ( DH * rho ) * A * p;
        vZ += v_factor * A * p;
        // velocity x: vX = vX + DT / ( DH * rho ) * B * p;
        vX += v_factor * B * p;
        // velocity y: vY = vY + DT / ( DH * rho ) * C * p;
        vY += v_factor * C * p;
        
        lama::Scalar znorm = vZ.l2Norm();
        lama::Scalar xnorm = vX.l2Norm();
        lama::Scalar ynorm = vY.l2Norm();
        
        // create new Vector(Pointer) with same configuration as vZ
        common::unique_ptr<lama::Vector> helpPtr( vZ.newVector() );
        // get Reference of VectorPointer
        lama::Vector& help = *helpPtr;
        
        // pressure update
        help =  DH_INV * D * vZ;
        help += DH_INV * E * vX;
        help += DH_INV * F * vY;
        p += p_factor * help; // p_factor is 'DT * M'
        
        lama::Scalar hnorm = help.l2Norm();
        lama::Scalar pnorm = p.l2Norm();
        
        // update seismogram and pressure with source terms
        // CAUTION: elementwise access by setVal and getVal cause performace issues executed on CUDA
        //          should be used rarely
        // TODO: can do this by index operator[] --> no need for DenseVector<>, can use Vector instead
        p.setValue( source_index, p.getValue( source_index ) + source.getValue( t ) );
        seismogram.setValue( t, p.getValue( seismogram_index ) );
        
        // TODO: plot snapshots of wave propagation ???
    }
}


