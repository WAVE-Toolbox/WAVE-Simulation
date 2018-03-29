#include <math.h>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "SourceSignal/SinThree.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(SinThreeTest, TestConstructor)
{
    int NT = 4;
    double DT = 1;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;

    lama::DenseVector<double> sampleResult;
    sampleResult.allocate(NT);

    //calculate sample result
    lama::DenseVector<double> zero(NT, 0.0);
    double temp;
    int time_index1, time_index2, i, count;
    time_index1 = floor(Tshift / DT);
    time_index2 = time_index1 + floor(1.0 / FC / DT);
    count = 0;
    for (i = time_index1; i <= time_index2; i++) {
        temp = count * DT * M_PI * FC;
        temp = sin(temp);
        temp = pow(temp, 3);
        zero.setValue(i, temp);
        count++;
    }
    sampleResult = lama::Scalar(AMP) * zero;

    //Testing
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    EXPECT_NO_THROW(Acquisition::SourceSignal::SinThree<double>(testResult1, NT, DT, FC, AMP, Tshift));

    lama::DenseVector<double> testResult2;
    testResult2.allocate(NT);
    Acquisition::SourceSignal::SinThree<double>(testResult2, NT, DT, FC, AMP, Tshift);

    EXPECT_EQ(sampleResult.getValue(3), testResult2.getValue(3));
}

TEST(SinThreeTest, TestAsserts)
{
    int NT = 4;
    double DT = 10;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);

    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinThree<double>(testResult1, -NT, DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinThree<double>(testResult1, NT, -DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::SinThree<double>(testResult1, NT, DT, -FC, AMP, Tshift));
}
