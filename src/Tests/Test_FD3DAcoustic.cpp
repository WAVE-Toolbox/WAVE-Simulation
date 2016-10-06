#include <scai/lama.hpp>
#include <iostream>


#include "../Configuration.hpp"

#define MASTER 0

#define HOST_PRINT( comm, msg )     \
{                                   \
int myRank = comm->getRank();   \
if ( myRank == MASTER )         \
{                               \
std::cout << msg;           \
}                               \
}


int main( int argc, char* argv[] )
{
    typedef double ValueType;
    
    if(argc!=2){
        std::cout<< "\n\nNo configuration file given!\n\n" << std::endl;
        return(2);
    }
    
    // read configuration parameter from file
    Configuration<ValueType> config(argv[1]);
    
    lama::DenseMatrix<ValueType> seismo_ref;
    lama::DenseMatrix<ValueType> seismo_syn;
    lama::DenseMatrix<ValueType> seismo_residual;
    scai::lama::Scalar L2_scalar=0.0;
    ValueType L2=0.0;
    
    seismo_ref.readFromFile(config.getSeismogramFilename());
    seismo_syn.readFromFile("ci/seismogram_ref.mtx");
    
    seismo_residual=(seismo_ref-seismo_syn);
    L2_scalar=seismo_residual.l2Norm();
    L2=L2_scalar.getValue<ValueType>();
    
    std::cout << "\n\nL2: " << L2 << std::endl;
    
    if(L2>0.01){
        std::cout << "Seismogram does not match reference solution.\n\n" << std::endl;
        return(1);
    }
    
    return 0;
}
