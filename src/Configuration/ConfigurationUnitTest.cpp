#include "Configuration.hpp"
#include "gtest/gtest.h"

using namespace KITGPI;

TEST(ConfigurationTest, readFromUnkownFile)
{
    ASSERT_ANY_THROW(Configuration::Configuration config("../src/Tests/Testfiles/configuration_100.txt"););
    Configuration::Configuration config("../src/Tests/Testfiles/configuration_1.txt");
    ASSERT_ANY_THROW(config.readFromFile("../src/Tests/Testfiles/configuration_100.txt"));
}

TEST(ConfigurationTest, getFunction)
{
    Configuration::Configuration config("../src/Tests/Testfiles/configuration_1.txt");
    ASSERT_EQ(1.2124445, config.get<double>("testvalue1"));
    ASSERT_EQ(100, config.get<int>("testvalue2"));
    ASSERT_EQ(-1.2124445, config.get<double>("testvalue3"));
    ASSERT_EQ(-100, config.get<int>("testvalue4"));
    ASSERT_EQ("test123", config.get<std::string>("testvalue5"));
    ASSERT_EQ("capiTAL", config.get<std::string>("TESTVALUE6"));
    ASSERT_EQ("capiTAL", config.get<std::string>("testvalue6"));
    ASSERT_NE("capital", config.get<std::string>("testvalue6"));
    ASSERT_ANY_THROW(config.get<std::string>("UnkownValue"));
    ASSERT_TRUE(config.get<bool>("testvalue7"));
    ASSERT_FALSE(config.get<bool>("testvalue8"));
    ASSERT_EQ("/file/path/test.mtx", config.get<std::string>("testvalue9"));
}

TEST(ConfigurationTest, readAdditionalConfigFromFile)
{
    Configuration::Configuration config("../src/Tests/Testfiles/configuration_1.txt");
    config.readFromFile("../src/Tests/Testfiles/configuration_2.txt");

    ASSERT_EQ(1.2124445, config.get<double>("testvalue1"));
    ASSERT_EQ(100, config.get<int>("testvalue2"));
    ASSERT_EQ(-1.2124445, config.get<double>("testvalue3"));
    ASSERT_EQ(-100, config.get<int>("testvalue4"));
    ASSERT_EQ("test123", config.get<std::string>("testvalue5"));
    ASSERT_EQ("capiTAL", config.get<std::string>("TESTVALUE6"));
    ASSERT_EQ("capiTAL", config.get<std::string>("testvalue6"));
    ASSERT_NE("capital", config.get<std::string>("testvalue6"));
    ASSERT_ANY_THROW(config.get<std::string>("UnkownValue"));
    ASSERT_TRUE(config.get<bool>("testvalue7"));
    ASSERT_FALSE(config.get<bool>("testvalue8"));
    ASSERT_EQ("/file/path/test.mtx", config.get<std::string>("testvalue9"));

    ASSERT_EQ("test/test/file.su", config.get<std::string>("additionalValue1"));
    ASSERT_EQ(123, config.get<int>("additionalValue2"));

    config.readFromFile("../src/Tests/Testfiles/configuration_2.txt", false);
    ASSERT_EQ(1.2124445, config.get<double>("testvalue1"));
    ASSERT_EQ(0.12345, config.get<double>("testvalue1"));
    ASSERT_EQ("test/test/file.su", config.get<std::string>("additionalValue1"));
    ASSERT_EQ(123, config.get<int>("additionalValue2"));

    config.readFromFile("../src/Tests/Testfiles/configuration_2.txt", true);
    ASSERT_NE(1.2124445, config.get<double>("testvalue1"));
    ASSERT_EQ(0.12345, config.get<double>("testvalue1"));
    ASSERT_EQ("test/test/file.su", config.get<std::string>("additionalValue1"));
    ASSERT_EQ(123, config.get<int>("additionalValue2"));
}

TEST(ConfigurationTest, print)
{
    Configuration::Configuration config("../src/Tests/Testfiles/configuration_2.txt");
    ASSERT_NO_THROW(config.print());
}
