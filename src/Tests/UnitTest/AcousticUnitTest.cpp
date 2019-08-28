//#define protected public
//#define private public

#include "Acoustic.hpp"
#include "ModelparameterFactory.hpp"
//#include "Factory.hpp"
#include "gtest/gtest.h"

using namespace KITGPI;
using namespace scai;

TEST(Acoustic, TestInit)
{
    Configuration::Configuration testConfig1("../src/Tests/Testfiles/configuration_4.txt");
    Configuration::Configuration testConfig2("../src/Tests/Testfiles/configuration_5.txt");
    hmemo::ContextPtr testCtx = hmemo::Context::getContextPtr();
    int numParameter = 4;
    dmemo::DistributionPtr no_dist_numParameter(new scai::dmemo::NoDistribution(numParameter));
    Modelparameter::Modelparameter<double>::ModelparameterPtr model(Modelparameter::Factory<double>::Create("acoustic"));
    Modelparameter::Acoustic<double> testObject;
    
    Acquisition::Coordinates<double> modelCoordinates(testConfig1);


    ASSERT_ANY_THROW(testObject.init(testConfig1, testCtx, no_dist_numParameter,modelCoordinates));
}
