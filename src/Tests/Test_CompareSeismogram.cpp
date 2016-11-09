#include <scai/lama.hpp>
#include <iostream>

#include "../Configuration/Configuration.hpp"
#include "../Common/HostPrint.hpp"

using namespace KITGPI;

int main( int argc, char* argv[] )
{
    typedef double ValueType;
    
    if(argc!=2){
        std::cout<< "\n\nNo configuration file given!\n\n" << std::endl;
        return(2);
    }
    
    // read configuration parameter from file
    Configuration::Configuration<ValueType> config(argv[1]);
    
    lama::DenseMatrix<ValueType> seismo_ref;
    lama::DenseMatrix<ValueType> seismo_syn;
    lama::DenseMatrix<ValueType> seismo_residual;
    scai::lama::Scalar L2_scalar=0.0;
    ValueType L2=0.0;
    
    std::string filenameRef=config.getSeismogramFilename();
    std::size_t pos=filenameRef.find(".ci.mtx");
    std::string filenameSyn = filenameRef.substr (0,pos) + ".ref.mtx";
    
    seismo_ref.readFromFile(filenameRef);
    seismo_syn.readFromFile(filenameSyn);
    
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