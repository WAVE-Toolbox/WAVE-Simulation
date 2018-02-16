#include <math.h>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "SourceSignal/FGaussian.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(FGaussianTest, TestConstructor)
{
    int NT = 4;
    double DT = 10;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;

    lama::DenseVector<double> sampleResult;
    sampleResult.allocate(NT);

    //calculate sample result
    lama::DenseVector<double> sampleT;
    sampleT.setRange(NT, double(0), DT);
    lama::DenseVector<double> sampleHelp(sampleT.size(), 1.2 / FC + Tshift);
    lama::DenseVector<double> sampleTau(sampleT - sampleHelp);
    sampleTau *= M_PI * FC;
    sampleHelp = -2.0 * sampleTau;
    sampleTau = -1.0 * sampleTau * sampleTau;
    sampleTau.exp();
    sampleResult = lama::Scalar(AMP) * sampleHelp * sampleTau;

    //Testing
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    EXPECT_NO_THROW(Acquisition::SourceSignal::FGaussian<double>(testResult1, NT, DT, FC, AMP, Tshift));

    lama::DenseVector<double> testResult2;
    testResult2.allocate(NT);
    Acquisition::SourceSignal::FGaussian<double>(testResult2, NT, DT, FC, AMP, Tshift);

    EXPECT_EQ(sampleResult.getValue(3), testResult2.getValue(3));
}

TEST(FGaussianTest, TestAsserts)
{
    int NT = 4;
    double DT = 10;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);

    EXPECT_ANY_THROW(Acquisition::SourceSignal::FGaussian<double>(testResult1, -NT, DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::FGaussian<double>(testResult1, NT, -DT, FC, AMP, Tshift));
    EXPECT_ANY_THROW(Acquisition::SourceSignal::FGaussian<double>(testResult1, NT, DT, -FC, AMP, Tshift));
}
