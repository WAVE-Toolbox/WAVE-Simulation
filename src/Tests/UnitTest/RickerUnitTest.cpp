#include <math.h>
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "SourceSignal/Ricker.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(RickerTest, TestConstructor)
{
    int NT = 4;
    double DT = 1;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;

    lama::DenseVector<double> sampleResult;
    sampleResult.allocate(NT);

    //calculate sample result
    auto t = lama::linearDenseVector<double>(NT, 0.0, DT);
    auto help = lama::fill<lama::DenseVector<double>>(t.size(), 1.5 / FC + Tshift);
    auto tau = lama::eval<lama::DenseVector<double>>(t - help);
    tau *= M_PI * FC;
    lama::DenseVector<double> one(sampleResult.size(), 1.0);
    help = tau * tau;
    tau = -1.0 * help;
    tau = lama::exp( tau );
    help = one - 2.0 * help;
    sampleResult = AMP * help * tau;

    //Testing
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);
    ASSERT_NO_THROW(Acquisition::SourceSignal::Ricker<double>(testResult1, NT, DT, FC, AMP, Tshift));

    lama::DenseVector<double> testResult2;
    testResult2.allocate(NT);
    Acquisition::SourceSignal::Ricker<double>(testResult2, NT, DT, FC, AMP, Tshift);

    ASSERT_EQ(sampleResult.getValue(3), testResult2.getValue(3));
}

TEST(RickerTest, TestAsserts)
{
    int NT = 4;
    double DT = 10;
    double FC = 1;
    double AMP = 1;
    double Tshift = 0;
    lama::DenseVector<double> testResult1;
    testResult1.allocate(NT);

    ASSERT_ANY_THROW(Acquisition::SourceSignal::Ricker<double>(testResult1, -NT, DT, FC, AMP, Tshift));
    ASSERT_ANY_THROW(Acquisition::SourceSignal::Ricker<double>(testResult1, NT, -DT, FC, AMP, Tshift));
    ASSERT_ANY_THROW(Acquisition::SourceSignal::Ricker<double>(testResult1, NT, DT, -FC, AMP, Tshift));
}
