#include <iostream>
#include <scai/common/Settings.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <string>
#include <vector>

#include "Configuration.hpp"
#include "HostPrint.hpp"

using namespace scai;
/*------------------------------
     TwoLayer - Model Creation
     This is a simple model creation tool which
     reads in the configuration file the model dimension and name
-------------------------------*/
int main(int argc, char *argv[])
{
    typedef double ValueType;

    // parameter
    ValueType vp1 = 3500, vs1 = 2000, rho1 = 2000, tauP1 = 0.1, tauS1 = 0.1;
    ValueType vp2 = 4550, vs2 = 2600, rho2 = 2600, tauP2 = 0.1, tauS2 = 0.1;

    // depth of second layer
    IndexType depth = 40;

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    // read configuration parameter from file
    KITGPI::Configuration::Configuration config(argv[1]);

    // estimate grid with parameters out of the configuration
    int NX = config.get<IndexType>("NX");
    int NY = config.get<IndexType>("NY");
    int NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NZ, NY, NX);

    // construct model vectors
    lama::GridVector<ValueType> vp(grid);
    lama::GridVector<ValueType> vs(grid);
    lama::GridVector<ValueType> rho(grid);
    lama::GridVector<ValueType> tauP(grid);
    lama::GridVector<ValueType> tauS(grid);

    //set velocities for first layer
    vp = vp1;
    vs = vs1;
    rho = rho1;
    tauP = tauP1;
    tauS = tauS1;

    //set velocities for second layer
    for (IndexType y = depth; y < NY; ++y) {

        vp(lama::Range(), y, lama::Range()) = vp2;
        vs(lama::Range(), y, lama::Range()) = vs2;
        rho(lama::Range(), y, lama::Range()) = rho2;
        tauP(lama::Range(), y, lama::Range()) = tauP2;
        tauS(lama::Range(), y, lama::Range()) = tauS2;
    }

    std::string type = config.get<std::string>("equationType");

    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");

    rho.writeToFile(filename + ".density.mtx");

    if (type.compare("sh") != 0) {
        vp.writeToFile(filename + ".vp.mtx");
    }

    if (type.compare("acoustic") != 0) {
        vs.writeToFile(filename + ".vs.mtx");
    }
    if (type.compare("visco") == 0) {
        tauP.writeToFile(filename + ".tauP.mtx");
        tauS.writeToFile(filename + ".tauS.mtx");
    }
    return 0;
}
