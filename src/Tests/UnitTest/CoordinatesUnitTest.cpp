#include <scai/lama.hpp>

#include "Coordinates.hpp"
#include <gtest/gtest.h>

using namespace scai;
using namespace KITGPI;

TEST(CoordinateTest, TestLocatedOnSurface)
{

    int NX = 5;
    int NY = 10;
    int NZ = 15;

    Acquisition::Coordinates test;

    int coord1 = 2;
    EXPECT_TRUE(test.locatedOnSurface(coord1, NX, NY, NZ));

    int coord2 = 52;
    EXPECT_TRUE(test.locatedOnSurface(coord2, NX, NY, NZ));

    int coord3 = 70;
    EXPECT_FALSE(test.locatedOnSurface(coord3, NX, NY, NZ));
}

TEST(CoordinateTest, TestIndex2coordinate)
{

    // Grid
    int NX = 5;
    int NY = 10;
    int NZ = 15;

    // Test coordinates
    int testX = 2;
    int testY = 2;
    int testZ = 2;

    // index of test coordinates
    int testCoord = 112;

    Acquisition::coordinate3D sampleSolution;
    Acquisition::coordinate3D result;

    Acquisition::Coordinates test;

    result = test.index2coordinate(testCoord, NX, NY, NZ);

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
    int NY = 10;
    int NZ = 15;

    // Test coordinates
    int testX = 4;
    int testY = 3;
    int testZ = 2;

    // index of test coordinates
    int sampleCordinate2index = 119;

    // Test first interface
    Acquisition::Coordinates test1;
    EXPECT_EQ(sampleCordinate2index, test1.coordinate2index(testX, testY, testZ, NX, NY, NZ));

    // Test second interface
    Acquisition::coordinate3D testCoord;
    testCoord.x = testX;
    testCoord.y = testY;
    testCoord.z = testZ;
    EXPECT_EQ(sampleCordinate2index, test1.coordinate2index(testCoord, NX, NY, NZ));

    // Test if interface throws if wrong input parameters are given
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, testZ, testX, NY, NZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, testZ, NX, testY, NZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, testZ, NX, NY, testZ));
    EXPECT_ANY_THROW(test1.coordinate2index(-testX, testY, testZ, NX, NY, NZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, -testY, testZ, NX, NY, NZ));
    EXPECT_ANY_THROW(test1.coordinate2index(testX, testY, -testZ, NX, NY, NZ));
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

    // Distance of these coordinates to edges
    Acquisition::coordinate3D solutionDistance;
    solutionDistance.x = 0;
    solutionDistance.y = 4;
    solutionDistance.z = 10;

    Acquisition::coordinate3D testCoord;
    testCoord.x = testX;
    testCoord.y = testY;
    testCoord.z = testZ;

    Acquisition::Coordinates test;
    Acquisition::coordinate3D result;
    result = test.edgeDistance(testCoord, NX, NY, NZ);

    EXPECT_EQ(solutionDistance.x, result.x);
    EXPECT_EQ(solutionDistance.y, result.y);
    EXPECT_EQ(solutionDistance.z, result.z);
}
