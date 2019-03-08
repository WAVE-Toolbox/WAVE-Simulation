#include "../Configuration/Configuration.hpp"
#include "Receivers.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

typedef float ValueType;

TEST(ReceiversTest, TestInit)
{

    Configuration::Configuration testConfig("../src/Tests/Testfiles/configuration_3.txt");

    // Create an object of the mapping (3D-1D) class Coordinates

    Acquisition::Coordinates<ValueType> testCoordinates(testConfig.get<IndexType>("NX"), testConfig.get<IndexType>("NY"), testConfig.get<IndexType>("NZ"), testConfig.get<ValueType>("DH"));

    hmemo::ContextPtr testCtx = hmemo::Context::getContextPtr();
    int numParameter = 4;
    dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(numParameter));
    Acquisition::Receivers<ValueType> testObject;

    ASSERT_ANY_THROW(testObject.init(testConfig, testCoordinates, testCtx, no_dist_numParameter));
}
