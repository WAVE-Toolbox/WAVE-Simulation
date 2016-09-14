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
 *  routine for calculation the input source (Ricker wavelet)
 *
 *  MATLAB:
 *  t=0:DT:(NT*DT-DT);
 *  tau=pi*FC*(t-1.5/FC);
 *  signal=AMP*(1-2*tau.^2).*exp(-tau.^2);
 */
template<typename ValueType>
void sourceFunction( lama::DenseVector<ValueType>& source, IndexType FC, IndexType AMP, dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "sourceFunction" )
    
    HOST_PRINT( comm, "Calculate source signal...\n" );
    
    // this is for tau[i] = pi * FC * ( source[i] - 1.5/FC );
    
    lama::DenseVector<ValueType> help( source.size(), 1.5 / FC );
    lama::DenseVector<ValueType> tau( source - help );
    tau *= M_PI * FC;
    
    // this is for source[i] = AMP * ( 1.0 - 2.0 * tau[i] * tau[i] * exp( -tau[i] * tau[i] ) );
    lama::DenseVector<ValueType> one( source.size(), 1.0 );
    help = tau * tau;
    tau = -1.0 * help;
    tau.exp();
    help = one - 2.0 * help;
    source = lama::Scalar(AMP) * help * tau;
}


