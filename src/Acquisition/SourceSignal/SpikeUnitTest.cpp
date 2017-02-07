#include <math.h>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "SinW.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(SpikeTest, TestConstructor)
{
    int NT=4;
    double DT=1;
    double FC=1;
    double AMP=1;
    double Tshift=0;
    
    lama::DenseVector<double> sampleResult;
    sampleResult.allocate(NT);
    
    //calculate sample result
    scai::lama::Scalar temp_spike;
    IndexType time_index;
    lama::DenseVector<double> help(NT, 0.0);
    temp_spike = 1.0;
    time_index = floor(Tshift / DT);
    help.setValue(time_index, temp_spike);
    sampleResult = lama::Scalar(AMP) * help;
    
    //Testing
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    EXPECT_NO_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, DT, FC, AMP, Tshift));
    
    lama::DenseVector<double> testResult2;
    testResult2.allocate(NT);
    Acquisition::SourceSignal::SinW<double>(testResult2, NT, DT, FC, AMP, Tshift);
    
    EXPECT_EQ(sampleResult.getValue(3),testResult2.getValue(3));
    
}

TEST(SpikeTest, TestAsserts)
{
    int NT=4;
    double DT=10;
    double FC=1;
    double AMP=1;
    double Tshift=0;
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    
    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, -NT, DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, -DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, DT, -FC, AMP, Tshift));
    
}
