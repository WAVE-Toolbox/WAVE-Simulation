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
                
        //! \brief Factory class.
        template<typename ValueType>
        class Factory
        {
            
        public:
            
            //! \brief Declare Wavefield pointer
            typedef typename Wavefields<ValueType>::WavefieldPtr  WavefieldPtr;

            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static WavefieldPtr Create( std::string dimension, std::string type ) {
                
                // transform to lower cases
                std::transform(type.begin(), type.end(), type.begin(), ::tolower);
                std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

                // Assert correctness of input values
                SCAI_ASSERT_ERROR(dimension.compare("2d")==0 ||dimension.compare("3d")==0, "Unkown dimension" );
                SCAI_ASSERT_ERROR(type.compare("acoustic")==0 || type.compare("elastic")==0 || type.compare("visco")==0, "Unkown type" );
                
                // 2D
                if( dimension.compare("2d") == 0 && type.compare("acoustic") == 0 ){
                    return WavefieldPtr(new FD2Dacoustic<ValueType>);
                }
                if( dimension.compare("2d") == 0 && type.compare("elastic") == 0 ){
                    return WavefieldPtr(new FD2Delastic<ValueType>);
                }
                if( dimension.compare("2d") == 0 && type.compare("visco") == 0 ){
                    return WavefieldPtr(new FD2Dvisco<ValueType>);
                }
                
                // 3D
                if( dimension.compare("3d") == 0 && type.compare("acoustic") == 0 ){
                    return WavefieldPtr(new FD3Dacoustic<ValueType>);
                }
                if( dimension.compare("3d") == 0 && type.compare("elastic") == 0 ){
                    return WavefieldPtr(new FD3Delastic<ValueType>);
                }
                if( dimension.compare("3d") == 0 && type.compare("visco") == 0 ){
                    return WavefieldPtr(new FD3Dvisco<ValueType>);
                }
                
                COMMON_THROWEXCEPTION("Reached end of factory without match");
                return nullptr;
            }
            
            
        };
    }
}
