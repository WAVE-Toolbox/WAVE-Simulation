#pragma once

#include <scai/lama/Scalar.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <string>
#include <algorithm> 
#include <cctype>
#include <scai/logging.hpp>

using namespace scai;

//! Namespace of the Geophysical Institute of the Karlsruhe Institute of Technology
namespace KITGPI {
    
    //! Configuration namespace
    namespace Configuration {
        
        //! Class for Configuration of the FD simulation
        /*!
         This class handels the configuration for the finite-difference simulation.
         */
        class Configuration
        {
        public:
            
            /*! \brief Default deconstructor
             */
            ~Configuration(){}
            
            explicit Configuration( std::string const& filename );
            
            void print() const;
            
            template<typename ValueType>
            ValueType get(std::string const& parameterName) const;
            
        private:
            
            std::unordered_map<std::string,std::string> configMap; ///< Map that is returned from config File
            
        };
    }
}


/*! \brief Constructor
 *
 \param filename of configuration file
 */
KITGPI::Configuration::Configuration::Configuration(std::string const& filename)
{
    // read all lines in file
    
    std::string line;
    std::ifstream input( filename.c_str() );
    while ( std::getline( input, line ) )
    {
        size_t lineEnd = line.size();
        std::string::size_type commentPos1 = line.find_first_of( "#", 0 );
        if ( std::string::npos != commentPos1 )
        {
            std::string::size_type commentPos2 = line.find_first_of( "#", commentPos1 );
            if ( std::string::npos != commentPos2 )
            {
                if( commentPos1 == 0 )
                {
                    continue;
                }
                lineEnd = commentPos1;
            }
        }
        
        std::string::size_type equalPos = line.find_first_of( "=", 0 );
        
        if ( std::string::npos != equalPos )
        {
            // tokenize it  name = val
            std::string name = line.substr( 0, equalPos );
            size_t len = lineEnd - ( equalPos + 1);
            std::string val  = line.substr( equalPos + 1, len);
            std::transform(name.begin(), name.end(), name.begin(), ::tolower);
            configMap.insert( std::pair<std::string,std::string>( name, val) );
        }
    }
    input.close();
}


/*! \brief Constructor
 *
 \param parameterName name of the parameter
 */
template<typename ValueType>
ValueType KITGPI::Configuration::Configuration::get( std::string const& parameterName) const
{
    ValueType temp;
    try {
        std::istringstream( configMap.at(parameterName) ) >> temp;
    }
    catch (...) {
        COMMON_THROWEXCEPTION("Parameter " << parameterName << ": Not found in Configuration file! " << std::endl)
    }
    return(temp);
}

/*! \brief Print configuration
 */
void KITGPI::Configuration::Configuration::print() const
{
    std::cout << "\t" << "Configuration: \n" << std::endl;
    for ( std::unordered_map<std::string,std::string>::const_iterator iter = configMap.begin(); iter != configMap.end(); ++iter )
        std::cout << "\t"<< iter->first << " = " << iter->second << std::endl;
    std::cout << std::endl;
}
