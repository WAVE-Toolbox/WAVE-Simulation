#include <iostream>
#include <scai/lama.hpp>

#include "HostPrint.hpp"
#include "Configuration.hpp"

int main(int argc, char *argv[])
{
    typedef double ValueType;

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);

    scai::lama::DenseMatrix<ValueType> seismo_ref;
    scai::lama::DenseMatrix<ValueType> seismo_syn;
    scai::lama::DenseMatrix<ValueType> seismo_residual;

    std::string filenameRef = config.get<std::string>("SeismogramFilename") + ".mtx";
    std::size_t pos = filenameRef.find(".ci.mtx");
    std::string filenameSyn = filenameRef.substr(0, pos) + ".ref.mtx";

    std::size_t found = filenameRef.find_last_of(".");
    std::string beforeEnding = filenameRef.substr(0, found);
    std::string afterEnding = filenameRef.substr(found);
    filenameRef = beforeEnding + ".p" + afterEnding;

    seismo_ref.readFromFile(filenameRef);
    seismo_syn.readFromFile(filenameSyn);

    seismo_residual = (seismo_ref - seismo_syn);
    auto L2 = seismo_residual.l2Norm();

    std::cout << "\n\nL2: " << L2 << std::endl;

    if (L2 > 0.01) {
        std::cout << "Seismogram does not match reference solution.\n\n"
                  << std::endl;
        return (1);
    } else {
        std::cout << "\n\n!!! Successful !!!\n\n"
                  << std::endl;
    }

    return 0;
}
