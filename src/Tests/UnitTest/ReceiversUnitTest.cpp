#define private public

#include "Receivers.hpp"
#include "../Configuration/Configuration.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

TEST(ReceiversTest, TestInit)
{
    
    Configuration::Configuration testConfig("../src/Tests/Testfiles/configuration_3.txt");
    
    scai::lama::DenseMatrix<double> acquisition_temp;
    acquisition_temp.readFromFile("../src/Tests/Testfiles/receiver.mtx");
    
    hmemo::ContextPtr testCtx = hmemo::Context::getContextPtr();
    int numParameter = 4;
    dmemo::DistributionPtr no_dist_numParameter( new scai::dmemo::NoDistribution ( numParameter ) );
    Acquisition::Receivers<double> testObject;
    
    
    ASSERT_ANY_THROW(testObject.init(acquisition_temp,testConfig, testCtx, no_dist_numParameter));
    
}


TEST(ReceiversTest, TestcheckRequiredNumParameter)
{
    Acquisition::Receivers<double> testObject;
    ASSERT_ANY_THROW(testObject.checkRequiredNumParameter(6));
    
}
