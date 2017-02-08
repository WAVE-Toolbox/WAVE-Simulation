#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "gtest/gtest.h"

using namespace scai;
using namespace KITGPI;

TEST(AcquisitionTest, TestParameters)
{

    ASSERT_EQ(0, (Acquisition::SeismogramType::P));
    ASSERT_EQ(1, (Acquisition::SeismogramType::VX));
    ASSERT_EQ(2, (Acquisition::SeismogramType::VY));
    ASSERT_EQ(3, (Acquisition::SeismogramType::VZ));
}

TEST(AcquisitionTest, TestMapping)
{

    ASSERT_EQ("p", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::P]));
    ASSERT_EQ("vx", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VX]));
    ASSERT_EQ("vy", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VY]));
    ASSERT_EQ("vz", (Acquisition::SeismogramTypeString[Acquisition::SeismogramType::VZ]));
}

TEST(AcquisitionTest, TestArrayLength)
{

    int nSeismogramType = sizeof(Acquisition::SeismogramType);
    EXPECT_EQ(nSeismogramType, static_cast<int>(Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE));
}
