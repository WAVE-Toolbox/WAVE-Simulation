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
    IndexType NX = config.get<IndexType>("NX");
    IndexType NY = config.get<IndexType>("NY");
    IndexType NZ = config.get<IndexType>("NZ");
    common::Grid3D grid(NY, NZ, NX);

    // construct model vectors and set velocities for first layer

    lama::GridVector<ValueType> vp(grid, vp1);
    lama::GridVector<ValueType> vs(grid, vs1);
    lama::GridVector<ValueType> rho(grid, rho1);
    lama::GridVector<ValueType> tauP(grid, tauP1);
    lama::GridVector<ValueType> tauS(grid, tauS1);

    //set velocities for second layer
    for (IndexType y = depth; y < NY; ++y) {

        vp(y, lama::Range(), lama::Range()) = vp2;
        vs(y, lama::Range(), lama::Range()) = vs2;
        rho(y, lama::Range(), lama::Range()) = rho2;
        tauP(y, lama::Range(), lama::Range()) = tauP2;
        tauS(y, lama::Range(), lama::Range()) = tauS2;
    }

    std::string type = config.get<std::string>("equationType");

    //write model to file specified in configuration
    std::string filename = config.get<std::string>("ModelFilename");

    std::string suffix = ".mtx"; // default output is matrix market

    IndexType fileFormat = config.get<IndexType>("FileFormat");
    if (fileFormat == 1) {

        rho.writeToFile(filename + ".density" + suffix);

        if (type.compare("sh") != 0) {
            vp.writeToFile(filename + ".vp" + suffix);
        }

        if (type.compare("acoustic") != 0) {
            vs.writeToFile(filename + ".vs" + suffix);
        }

        if (type.compare("visco") == 0) {
            tauP.writeToFile(filename + ".tauP" + suffix);
            tauS.writeToFile(filename + ".tauS" + suffix);
        }

    } else {
        if (fileFormat == 2)
            suffix = ".lmf";
        if (fileFormat == 3)
            suffix = ".frv";

        rho.writeToFile(filename + ".density" + suffix, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);

        if (type.compare("sh") != 0) {
            vp.writeToFile(filename + ".vp" + suffix, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
        }

        if (type.compare("acoustic") != 0) {
            vs.writeToFile(filename + ".vs" + suffix, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
        }
        if (type.compare("visco") == 0) {
            tauP.writeToFile(filename + ".tauP" + suffix, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
            tauS.writeToFile(filename + ".tauS" + suffix, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
        }
    }

    return 0;
}
