#include <math.h>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "SourceSignal/SinW.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(SpikeTest, TestConstructor)
{
    int NT = 4;
    double DT = 1;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;

    lama::DenseVector<double> sampleResult;
    sampleResult.allocate(NT);

    //calculate sample result
    IndexType time_index;
    lama::DenseVector<double> help(NT, 0.0);
    double temp_spike = 1.0;
    time_index = floor(Tshift / DT);
    help.setValue(time_index, temp_spike);
    sampleResult = AMP * help;

    //Testing
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    ASSERT_NO_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, DT, FC, AMP, Tshift));

    lama::DenseVector<double> testResult2;
    testResult2.allocate(NT);
    Acquisition::SourceSignal::SinW<double>(testResult2, NT, DT, FC, AMP, Tshift);

    ASSERT_EQ(sampleResult.getValue(3), testResult2.getValue(3));
}

TEST(SpikeTest, TestAsserts)
{
    int NT = 4;
    double DT = 10;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);

    ASSERT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, -NT, DT, FC, AMP, Tshift));
    ASSERT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, -DT, FC, AMP, Tshift));
    ASSERT_ANY_THROW(Acquisition::SourceSignal::SinW<double>(testResult1, NT, DT, -FC, AMP, Tshift));
}
