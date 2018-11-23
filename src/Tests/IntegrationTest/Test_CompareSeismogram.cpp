#include <iostream>
#include <scai/common/Settings.hpp>
#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "Configuration.hpp"
#include "HostPrint.hpp"
#include "Receivers.hpp"

using namespace KITGPI;
using namespace scai;

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

    /* --------------------------------------- */
    /* Context and Distribution                */
    /* --------------------------------------- */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT
    dmemo::DistributionPtr dist(new scai::dmemo::NoDistribution(config.get<IndexType>("NX") * config.get<IndexType>("NY") * config.get<IndexType>("NZ")));

    KITGPI::Acquisition::Seismogram<ValueType> seismogramTest;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramRef;
    KITGPI::Acquisition::Seismogram<ValueType> seismogramDiff;

    Acquisition::Receivers<ValueType> receiversTest;
    receiversTest.init(config, ctx, dist);
    Acquisition::Receivers<ValueType> receiversRef;
    receiversRef.init(config, ctx, dist);

    std::string filenameTest = config.get<std::string>("SeismogramFilename") + ".mtx";
    std::string filenameRef = config.get<std::string>("SeismogramFilename") + ".mtx";
    std::size_t pos = filenameRef.find(".ci.mtx");
    filenameRef = filenameRef.substr(0, pos) + ".ref.mtx";

    receiversTest.getSeismogramHandler().read(config,filenameTest, 1);
    receiversRef.getSeismogramHandler().read(config,filenameRef, 1);

    ValueType misfit = 0, misfitSum = 0;

    for (int i = 0; i < KITGPI::Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismogramTest = receiversTest.getSeismogramHandler().getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));
        seismogramRef = receiversRef.getSeismogramHandler().getSeismogram(static_cast<KITGPI::Acquisition::SeismogramType>(i));

        seismogramDiff = seismogramTest - seismogramRef;
        misfit = seismogramDiff.getData().l2Norm();
        misfitSum += misfit;
    }

    std::cout << "\n\nL2: " << misfitSum << std::endl;

    if (misfitSum > 0.01) {
        std::cout << "Seismogram does not match reference solution.\n\n"
                  << std::endl;
        return (1);
    } else {
        std::cout << "\n\n!!! Successful !!!\n\n"
                  << std::endl;
    }

    return 0;
}
