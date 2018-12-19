#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include "BoundaryCondition/ABS2D.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

bool verbose; // global variable definition (needed in only one unit test)

TEST(ABS2DTest, TestApplyThrows)
{
    int N = 10;
    double testValue = 123.0;
    lama::DenseVector<double> testVector;
    testVector.allocate(N);
    testVector=testValue;
    
    ForwardSolver::BoundaryCondition::ABS2D<double> test;
    
    ASSERT_ANY_THROW(test.apply(testVector, testVector, testVector));
    ASSERT_ANY_THROW(test.apply(testVector, testVector, testVector, testVector, testVector));
    
}
