#pragma once

#include <cctype>
#include <fstream>
#include <iostream>
#include <scai/lama.hpp>
#include <sstream>
#include <string>
#include <unordered_map>

using namespace scai;

//! Namespace of the Geophysical Institute of the Karlsruhe Institute of Technology
namespace KITGPI
{

    //! Configuration namespace
    namespace Configuration
    {

        //! Class to handle the configuration.
        /*!
         * This class handels all aspects regarding reading and getting the configuration from a text file.
         * See the readFromFile() member methhod for information how the text file is handled.
         *
         */
        class Configuration
        {
          public:
            //! Default destructor
            ~Configuration(){};

            explicit Configuration(std::string const &filename);

            void readFromFile(std::string const &filename, bool overwrite);

            void print() const;

            template <typename ReturnType>
            ReturnType get(std::string const &parameterName) const;

          private:
            void add2map(std::string const &KEY, std::string const &VALUE, bool overwrite);

            std::unordered_map<std::string, std::string> configMap; ///< Map that is used for internal handling of the `KEY=VALUE` pairs
        };
    }
}

/*! \brief Add `KEY=VALUE` pairs to the configuration
 *
 * This function adds `KEY=VALUE` pairs to the configuration class.
 *
 * If a `KEY` is already set in the configuration, the `VALUE` of this `KEY` will not be overwritten by the new `VALUE`.
 * However, this behavior can be changed by the bool `overwrite`. If `overwrite` is set to `true` existing `VALUE`s will be overwritten.
 * The default value for `overwrite` is `false`.
 *
 \param KEY which correspondents to the VALUE
 \param VALUE for the KEY to add
 \param overwrite Bool which indicates if `VALUE`s for existing `KEY`s will be overwritten by the new `VALUE`
 */
void KITGPI::Configuration::Configuration::add2map(std::string const &KEY, std::string const &VALUE, bool overwrite = false)
{
    if (configMap.count(KEY) == 0) {
        configMap.insert(std::pair<std::string, std::string>(KEY, VALUE));
    } else {
        if (overwrite) {
            configMap.erase(KEY);
            configMap.insert(std::pair<std::string, std::string>(KEY, VALUE));
        }
    }
}

/*! \brief Read a configuration text file and adds its content to the configuration
 *
 * The text file is assumed to have the form:\n
 * `KEY=VALUE`\n
 * The `KEY` and `VALUE` are internaly handled as `std::string`. 
 * However, the `VALUE` will be casted to the requested data type by the get() function.
 * See the get() function for more information regarding getting a `VALUE` for `KEY` and for information on casting.
 *
 * The configuration file can also contain comments, which should be indicated by "#". Two types of comments are supported:
 * 1. Whole lines of comments e.g:\n
 *  `# This whole line is a comment`
 * 2. Comments after a KEY-VALUE pair, e.g.:\n
 *  `KEY=VALUE # comment`
 *
 \param filename of the configuration file to read in
 \param overwrite Bool which indicates if existing entries will be overriden or not
 */
void KITGPI::Configuration::Configuration::readFromFile(std::string const &filename, bool overwrite = false)
{
    std::string line;
    std::ifstream input(filename.c_str());
    if (input.good() != true) {
        COMMON_THROWEXCEPTION("Configuration file " << filename << " was not found " << std::endl)
    }
    while (std::getline(input, line)) {
        size_t lineEnd = line.size();
        std::string::size_type commentPos1 = line.find_first_of("#", 0);
        if (std::string::npos != commentPos1) {
            std::string::size_type commentPos2 = line.find_first_of("#", commentPos1);
            if (std::string::npos != commentPos2) {
                if (commentPos1 == 0) {
                    continue;
                }
                lineEnd = commentPos1;
            }
        }

        std::string::size_type equalPos = line.find_first_of("=", 0);

        if (std::string::npos != equalPos) {
            // tokenize it  name = val
            std::string name = line.substr(0, equalPos);
            size_t len = lineEnd - (equalPos + 1);
            std::string val = line.substr(equalPos + 1, len);
            std::transform(name.begin(), name.end(), name.begin(), ::tolower);

            add2map(name, val, overwrite);
        }
    }
    input.close();
}

/*! \brief Constructor which reads in the configuration from file
 *
 * The constructor reads in the configuration from a text file.\n
 *
 * See the readFromFile() member function regarding information how the text file is parsed.
 *
 \param filename of the configuration file to read in
 */
KITGPI::Configuration::Configuration::Configuration(std::string const &filename)
{
    readFromFile(filename, true);
}

/*! \brief Get the value of a parameter from the read-in configuration
 *
 * This method returns the value of a specific parameter from the configuration.
 * If the specified value is not found in the configuration, e.g. it was not included in the
 * configuration file, the function will throw an exception. The exception message will contain the missing parameter. \n
 * This getter function is tested with the following data types: int (IndexType), float, double (ValueType) and std::string.
 *
 * **Known limitations:**\n
 * The conversation to bool fails if strings are used as `VALUE`.\n
 * This will fail: `Test1=true` -> `config.get<bool>("Test1") == true;`\n
 * However, this will work:`Test1=1` -> `config.get<bool>("Test1") == true;`
 *
 * **Usage:**\n
 * `int test = config.get<int>("test"); // For an int` \n
 * `std::string filename = config.get<std::string>("filename"); // for a std::string` \n
 * `ReturnType VALUE = config.get<ReturnType>("KEY"); // general`
 \param parameterName name of the wanted parameter
 \return The value of "parameterName"
 */
template <typename ReturnType>
ReturnType KITGPI::Configuration::Configuration::get(std::string const &parameterName) const
{
    std::string tempName = parameterName;
    ReturnType temp;
    try {
        std::transform(tempName.begin(), tempName.end(), tempName.begin(), ::tolower);
        std::istringstream(configMap.at(tempName)) >> temp;
    } catch (...) {
        COMMON_THROWEXCEPTION("Parameter " << parameterName << ": Not found in Configuration file! " << std::endl)
    }
    return (temp);
}

/*! \brief Print configuration to stdout
 *
 * This function prints all `KEY=VALUE` pairs contained in the read-in configuration to stdout.
 */
void KITGPI::Configuration::Configuration::print() const
{
    std::cout << "\t"
              << "Configuration: \n"
              << std::endl;
    for (std::unordered_map<std::string, std::string>::const_iterator iter = configMap.begin(); iter != configMap.end(); ++iter)
        std::cout << "\t" << iter->first << " = " << iter->second << std::endl;
    std::cout << std::endl;
}
