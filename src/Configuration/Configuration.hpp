#pragma once

#include <cctype>
#include <fstream>
#include <iostream>
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

            void readFromFile(std::string const &filename, bool overwrite = false);

            void print() const;

            template <typename ReturnType>
            ReturnType get(std::string const &parameterName) const;
            void add2config(std::string const &KEY, std::string const &VALUE);
            
          private:
            void add2map(std::string const &KEY, std::string const &VALUE, bool overwrite = false);

            std::unordered_map<std::string, std::string> configMap; ///< Map that is used for internal handling of the `KEY=VALUE` pairs
        };
    }
}
