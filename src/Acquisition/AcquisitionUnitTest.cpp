#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(AcquisitionTest, Parameters)
{

    EXPECT_EQ(0, (Acquisition::SeismogramType::P));
    EXPECT_EQ(1, (Acquisition::SeismogramType::VX));
    EXPECT_EQ(2, (Acquisition::SeismogramType::VY));
    EXPECT_EQ(3, (Acquisition::SeismogramType::VZ));
}

TEST(AcquisitionTest, Mapping)
{

    EXPECT_EQ("p", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::P]));
    EXPECT_EQ("vx", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VX]));
    EXPECT_EQ("vy", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VY]));
    EXPECT_EQ("vz", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VZ]));
}

TEST(AcquisitionTest, ArrayLength)
{

    int nSeismogramType = sizeof(Acquisition::SeismogramType);
    EXPECT_EQ(nSeismogramType, static_cast<int>(Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE));
}
