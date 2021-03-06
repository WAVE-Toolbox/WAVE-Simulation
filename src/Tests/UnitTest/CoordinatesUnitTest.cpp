#include <scai/lama.hpp>

#include "Coordinates.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

typedef float ValueType;

TEST(CoordinateTest, TestLocatedOnSurface)
{
    IndexType NX = 5;
    IndexType NY = 10;
    IndexType NZ = 15;
    ValueType DH = 1.0;

    Acquisition::Coordinates<ValueType> test(NX, NY, NZ, DH);

    IndexType coord1 = 2;
    EXPECT_TRUE(test.locatedOnSurface(coord1));
    IndexType coord2 = 52;
    EXPECT_TRUE(test.locatedOnSurface(coord2));
    IndexType coord3 = 80;
    EXPECT_FALSE(test.locatedOnSurface(coord3));
}

TEST(CoordinateTest, TestIndex2coordinate)
{

    // Grid
    IndexType NX = 5;
    IndexType NY = 15;
    IndexType NZ = 10;

    // Test coordinates
    IndexType testX = 2;
    IndexType testY = 2;
    IndexType testZ = 2;
    double DH = 1.0;

    // index of test coordinates
    IndexType testIndex= 112;

    Acquisition::coordinate3D sampleSolution;
    Acquisition::coordinate3D result;

    Acquisition::Coordinates<ValueType> test(NX, NY, NZ, DH);

    result = test.index2coordinate(testIndex);

    sampleSolution.x = testX;
    sampleSolution.y = testY;
    sampleSolution.z = testZ;

    EXPECT_EQ(sampleSolution.x, result.x);
    EXPECT_EQ(sampleSolution.y, result.y);
    EXPECT_EQ(sampleSolution.z, result.z);
}

TEST(CoordinateTest, TestCoordinate2index)
{

    // Grid
    int NX = 5;
    int NY = 15;
    int NZ = 10;

    // Test coordinates
    int testX = 4;
    int testY = 2;
    int testZ = 3;
    double DH = 1.0;

    // index of test coordinates
    int sampleCordinate2index = 119;

    // Test first interface
    Acquisition::Coordinates<ValueType> test1(NX, NY, NZ, DH);
    EXPECT_EQ(sampleCordinate2index, test1.coordinate2index(testX, testY, testZ));

    // Test second interface
    Acquisition::coordinate3D testCoord;
    testCoord.x = testX;
    testCoord.y = testY;
    testCoord.z = testZ;
    EXPECT_EQ(sampleCordinate2index, test1.coordinate2index(testCoord));

    // Test if interface throws if wrong input parameters are given
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, 10 * testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, 8 * testY, testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(2 * testX, testY, testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(-testX, testY, testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, -testY, testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, -testZ));
}

TEST(CoordinateTest, TestEstimateDistanceToEdges3D)
{

    // Grid
    int NX = 5;
    int NY = 10;
    int NZ = 15;

    // Test coordinates
    int testX = 4;
    int testY = 4;
    int testZ = 4;
    double DH = 1.0;

    // Distance of these coordinates to edges
    Acquisition::coordinate3D solutionDistance;
    solutionDistance.x = 0;
    solutionDistance.y = 4;
    solutionDistance.z = 4;

    Acquisition::coordinate3D testCoord;
    testCoord.x = testX;
    testCoord.y = testY;
    testCoord.z = testZ;

    Acquisition::Coordinates<ValueType> test(NX, NY, NZ, DH);
    Acquisition::coordinate3D result;
    result = test.edgeDistance(testCoord);

    EXPECT_EQ(solutionDistance.x, result.x);
    EXPECT_EQ(solutionDistance.y, result.y);
    EXPECT_EQ(solutionDistance.z, result.z);
}
