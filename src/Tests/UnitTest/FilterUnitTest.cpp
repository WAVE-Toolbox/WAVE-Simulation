#include "../../Filter/Filter.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

typedef double ValueType;

TEST(FilterTest, TestLowpass)
{
    lama::DenseMatrix<ValueType> signal;
    signal.readFromFile("../src/Tests/Testfiles/filterTest_signal.mtx");
    
    lama::DenseMatrix<ValueType> signalFiltRef;
    signalFiltRef.readFromFile("../src/Tests/Testfiles/filterTest_signalFiltRef.mtx");
    
    Filter::Filter<ValueType> freqFilter;
    freqFilter.init(1e-3, 500);
    freqFilter.calc("butterworth", "lp", 4, 15.0);
    freqFilter.apply(signal);
    
    signal -= signalFiltRef;
    EXPECT_LT(signal.l2Norm(), 0.005);
}