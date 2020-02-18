
#pragma once

#include <cctype>
#include <fstream>
#include <iostream>
#include <list>
#include <scai/lama.hpp>
#include <sstream>
#include <string>
#include <unordered_map>

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
            //! Default constructor
            Configuration(){};

            explicit Configuration(std::string const &filename);

            void init(std::string const &filename);
            
            void readFromFile(std::string const &filename, bool overwrite = false);

            void print() const;

            template <typename ReturnType>
            ReturnType get(std::string const &parameterName) const;

            template <typename InputType>
            void add2config(std::string const &KEY, InputType const &VALUE, bool overwrite = false);

          private:
            void add2map(std::string const &KEY, std::string const &VALUE, bool overwrite = false);

            std::unordered_map<std::string, std::string> configMap; ///< Map that is used for internal handling of the `KEY=VALUE` pairs
            std::list<std::string> insertionOrder;                  ///< list of insertion order
        };

        /*! \brief Add `KEY=VALUE` pairs to the configuration
	*
	* This function adds `KEY=VALUE` pairs to the configuration class.
	*
	\param KEY which correspondents to the VALUE
	\param VALUE for the KEY to add
	*/
        template <typename InputType>
        void Configuration::add2config(std::string const &KEY, InputType const &VALUE, bool overwrite)
        {
            std::string tempName = KEY;
            std::transform(tempName.begin(), tempName.end(), tempName.begin(), ::tolower);

            std::ostringstream sstream;
            sstream << VALUE;
            std::string tempValue = sstream.str();

            add2map(tempName, tempValue, overwrite);
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
                std::istringstream input(configMap.at(tempName));
                input >> temp;
            } catch (...) {
                COMMON_THROWEXCEPTION("Parameter " << parameterName << ": Not found in Configuration file! " << std::endl)
            }
            return (temp);
        }
    }
}
