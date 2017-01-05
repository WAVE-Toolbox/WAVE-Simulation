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

SCAI_LOG_DEF_TEMPLATE_LOGGER(template<typename ValueType>, KITGPI::Configuration::Configuration<ValueType>::logger, "Configuration" )

//! Namespace of the Geophysical Institute of the Karlsruhe Institute of Technology
namespace KITGPI {
    
    //! Configuration namespace
    namespace Configuration {
        
        //! Class for Configuration of the FD simulation
        /*!
         This class handels the configuration for the finite-difference simulation.
         */
        template<typename ValueType>
        class Configuration
        {
        public:
            
            /*! \brief Default deconstructor
             */
            ~Configuration(){}
            
            explicit Configuration( std::string const& filename );
            
            void print() const;
            
            bool checkConfigPlausibility();
            
            std::unordered_map<std::string,std::string>& getMap() const { return configMap; } ///< Return configuration map
            
            IndexType getIndex(std::string const& parameterName) const;
            ValueType getValue(std::string const& parameterName) const;
            std::string getString(std::string const& parameterName) const;
            
        private:
            SCAI_LOG_DECL_STATIC_LOGGER(logger);
            
            std::unordered_map<std::string,std::string> configMap; ///< Map that is returned from config File
            
        };
    }
}


/*! \brief Constructor
 *
 \param filename of configuration file
 */
template<typename ValueType>
KITGPI::Configuration::Configuration<ValueType>::Configuration(std::string const& filename)
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
    SCAI_LOG_DEBUG(logger, "Map has been created");
    input.close();
}


/*! \brief Constructor
 *
 \param parameterName name of the parameter
 */
template<typename ValueType>
IndexType KITGPI::Configuration::Configuration<ValueType>::getIndex( std::string const& parameterName) const
{
    IndexType temp(0);
    try {
        std::string parameterNameTemp = parameterName;
        std::transform(parameterNameTemp.begin(), parameterNameTemp.end(), parameterNameTemp.begin(), ::tolower);
        std::istringstream( configMap.at(parameterNameTemp) ) >> temp;
    }
    catch (...) {
        COMMON_THROWEXCEPTION( "Parameter " << parameterName << ": Not found in configuration file! " << std::endl)
    }
    return static_cast<IndexType>(temp);
}


/*! \brief Constructor
 *
 \param parameterName name of the parameter
 */
template<typename ValueType>
ValueType KITGPI::Configuration::Configuration<ValueType>::getValue( std::string const& parameterName) const
{
    ValueType temp(0);
    try {
        std::string parameterNameTemp = parameterName;
        std::transform(parameterNameTemp.begin(), parameterNameTemp.end(), parameterNameTemp.begin(), ::tolower);
        std::istringstream( configMap.at(parameterNameTemp) ) >> temp;
    }
    catch (...) {
        COMMON_THROWEXCEPTION("Parameter " << parameterName << ": Not found in configuration file! " << std::endl)
    }
    return static_cast<ValueType>(temp);
}


/*! \brief Constructor
 *
 \param parameterName name of the parameter
 */
template<typename ValueType>
std::string KITGPI::Configuration::Configuration<ValueType>::getString( std::string const& parameterName) const
{
    std::string temp;
    try {
        std::string parameterNameTemp = parameterName;
        std::transform(parameterNameTemp.begin(), parameterNameTemp.end(), parameterNameTemp.begin(), ::tolower);
        temp = std::istringstream( configMap.at(parameterNameTemp) ).str();
    }
    catch (...) {
    COMMON_THROWEXCEPTION("String " << parameterName << ": Not found in configuration file! " )
    }
    return temp;
    
}


/*! \brief Print configuration
 */
template<typename ValueType>
void KITGPI::Configuration::Configuration<ValueType>::print() const
{
    std::cout << "\t" << "Configuration: \n" << std::endl;
    for ( std::unordered_map<std::string,std::string>::const_iterator iter = configMap.begin(); iter != configMap.end(); ++iter )
    std::cout << "\t"<< iter->first << " = " << iter->second << std::endl;
    std::cout << std::endl;
}


/*! \brief Check plausibility of the configuration file
 */
template<typename ValueType>
bool KITGPI::Configuration::Configuration<ValueType>::checkConfigPlausibility()
{
    if ((getIndex("spatialFDorder") % 2 != 0) || (getIndex("spatialFDorder") < 2) || (getIndex("spatialFDorder") > 12))
    {
        return false;
    }
    return true;	// check is positive
}
