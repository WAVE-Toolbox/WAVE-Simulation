#pragma once

#include "Wavefields2Dacoustic.hpp"
#include "Wavefields2Delastic.hpp"
#include "Wavefields2Dvisco.hpp"
#include "Wavefields3Dacoustic.hpp"
#include "Wavefields3Delastic.hpp"
#include "Wavefields3Dvisco.hpp"
#include "Wavefields.hpp"

namespace KITGPI {
    
    //! \brief Wavefields namespace
    namespace Wavefields{
                

        template<typename ValueType>
        class Factory
        {
            
        public:
            
            typedef typename Wavefields<ValueType>::WavefieldPtr  WavefieldPtr;

            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            static WavefieldPtr Create( std::string const& dimension, std::string const& type ) {
                
                // Assert correctness of input values
                SCAI_ASSERT_ERROR(dimension.compare("2D")==0 ||dimension.compare("3D")==0, "Unkown dimension" );
                SCAI_ASSERT_ERROR(type.compare("acoustic")==0 || type.compare("elastic")==0 || type.compare("visco")==0, "Unkown type" );
                
                // 2D
                if( dimension.compare("2D") == 0 && type.compare("acoustic") == 0 ){
                    return WavefieldPtr(new FD2Dacoustic<ValueType>);
                }
                if( dimension.compare("2D") == 0 && type.compare("elastic") == 0 ){
                    return WavefieldPtr(new FD2Delastic<ValueType>);
                }
                if( dimension.compare("2D") == 0 && type.compare("visco") == 0 ){
                    return WavefieldPtr(new FD2Dvisco<ValueType>);
                }
                
                // 3D
                if( dimension.compare("3D") == 0 && type.compare("acoustic") == 0 ){
                    return WavefieldPtr(new FD3Dacoustic<ValueType>);
                }
                if( dimension.compare("3D") == 0 && type.compare("elastic") == 0 ){
                    return WavefieldPtr(new FD3Delastic<ValueType>);
                }
                if( dimension.compare("3D") == 0 && type.compare("visco") == 0 ){
                    return WavefieldPtr(new FD3Dvisco<ValueType>);
                }
                
                COMMON_THROWEXCEPTION("Reached end of factory without match");
                return nullptr;
            }
            
            
        };
    }
}
