#include "AcquisitionGeometry.hpp"
#include "Receivers.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

TEST(AcquisitionGeometryTest, TestThrowReadAcquisitionFromFile)
{
    int NX = 10;
    int NY = 20;
    int NZ = 30;
    
    std::string filename = "testReceiverFilename";
    hmemo::ContextPtr testCtx = hmemo::Context::getContextPtr();
    int numParameter = 4;
    dmemo::DistributionPtr no_dist_numParameter( new scai::dmemo::NoDistribution ( numParameter ) );
    
    Acquisition::Receivers<double> testObject;
    
    ASSERT_ANY_THROW(testObject.readAcquisitionFromFile(filename,-NX,NY,NZ,no_dist_numParameter,testCtx));
    ASSERT_ANY_THROW(testObject.readAcquisitionFromFile(filename,NX,-NY,NZ,no_dist_numParameter,testCtx));
    ASSERT_ANY_THROW(testObject.readAcquisitionFromFile(filename,NX,NY,-NZ,no_dist_numParameter,testCtx));
    
}



