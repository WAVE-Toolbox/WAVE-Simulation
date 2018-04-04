#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "BoundaryCondition/ABS3D.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(ABS3DTest, TestApplyThrows)
{
    int N = 10;
    double testValue = 123.0;
    lama::DenseVector<double> testVector;
    testVector.allocate(N);
    testVector=testValue;
    
    ForwardSolver::BoundaryCondition::ABS3D<double> test;
    
    ASSERT_ANY_THROW(test.apply(testVector, testVector, testVector, testVector));
    ASSERT_ANY_THROW(test.apply(testVector, testVector, testVector, testVector, testVector, testVector, testVector, testVector, testVector));
    
}
