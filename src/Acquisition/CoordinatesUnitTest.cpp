#include <scai/lama.hpp>
#include "gtest/gtest.h"
#include "Coordinates.hpp"

using namespace scai;
using namespace KITGPI;

TEST(CoordinateTest, locatedOnSurface){
    
    int NX = 5;
    int NY = 10;
    int NZ = 15;
    
    Acquisition::Coordinates<double> test;
    
    int coord1 = 2;
    EXPECT_TRUE(test.locatedOnSurface(coord1,NX,NY,NZ));
    
    int coord2 = 52;
    EXPECT_TRUE(test.locatedOnSurface(coord2,NX,NY,NZ));
    
    int coord3 = 70;
    EXPECT_FALSE(test.locatedOnSurface(coord3,NX,NY,NZ));
}


//TEST(CoordinateTest, index2coordinate){
//
//    int NX = 5;
//    int NY = 10;
//    int NZ = 15;
//    int testX = 5;
//    int testY = 4;
//    int testZ = 3;
//    
//    int testCoord = ( ( testX ) + ( testY ) * NX + ( testZ ) * NX * NY );
//    
//    Acquisition::coordinate3D sampleSolution;
//    Acquisition::Coordinates<double> test;
//    Acquisition::coordinate3D result;
//    result = test.index2coordinate(testCoord,NX,NY,NZ);
//    
//    sampleSolution.x = testX;
//    sampleSolution.y = testY;
//    sampleSolution.z = testZ;
//    
//    EXPECT_EQ(sampleSolution.x, result.x);
//    EXPECT_EQ(sampleSolution.y, result.y);
//    EXPECT_EQ(sampleSolution.z, result.z);
//    
//}
//
//
//TEST(CoordinateTest, coordinate2index){
//    
//    int NX = 5;
//    int NY = 10;
//    int NZ = 15;
//    int testX = 5;
//    int testY = 4;
//    int testZ = 3;
//    
//    Acquisition::Coordinates<double> test;
//    int sampleCordinate2index = ( ( testX ) + ( testY ) * NX + ( testZ ) * NX * NY );
//    
//    EXPECT_EQ(sampleCordinate2index, test.coordinate2index(testX,testY,testZ,NX,NY,NZ));
//    EXPECT_FALSE(test.coordinate2index(testX,testY,testZ,testX,NY,NZ));
//    EXPECT_FALSE(test.coordinate2index(testX,testY,testZ,NX,testY,NZ));
//    EXPECT_FALSE(test.coordinate2index(testX,testY,testZ,NX,NY,testZ));
//    EXPECT_FALSE(test.coordinate2index(-testX,testY,testZ,NX,NY,NZ));
//    EXPECT_FALSE(test.coordinate2index(testX,-testY,testZ,NX,NY,NZ));
//    EXPECT_FALSE(test.coordinate2index(testX,testY,-testZ,NX,NY,NZ));
//    
//    Acquisition::coordinate3D testCoord;
//    testCoord.x = testX;
//    testCoord.y = testY;
//    testCoord.z = testZ;
//    
//    Acquisition::Coordinates<double> test2;
//    
//    EXPECT_EQ(sampleCordinate2index, test2.coordinate2index(testCoord,NX,NY,NZ));
//    
//}
//
//
//TEST(CoordinateTest, estimateDistanceToEdges3D){
//    
//    int NX = 5;
//    int NY = 10;
//    int NZ = 15;
//    int testX = 4;
//    int testY = 4;
//    int testZ = 4;
//    
//    Acquisition::coordinate3D solutionDistance;
//    
//    solutionDistance.x=!((NX-testX)<(testX-1))?(testX-1):(NX-testX);
//    solutionDistance.y=!((NY-testY)<(testY-1))?(testY-1):(NY-testY);
//    solutionDistance.z=!((NX-testZ)<(testZ-1))?(testZ-1):(NZ-testZ);
//    
//    Acquisition::coordinate3D testCoord;
//    testCoord.x = testX;
//    testCoord.y = testY;
//    testCoord.z = testZ;
//    
//    Acquisition::Coordinates<double> test;
//    Acquisition::coordinate3D result;
//    result = test.edgeDistance(testCoord,NX,NY,NZ);
//    
//    EXPECT_EQ(solutionDistance.x, result.x);
//    EXPECT_EQ(solutionDistance.y, result.y);
//    EXPECT_EQ(solutionDistance.z, result.z);
//    
//}
