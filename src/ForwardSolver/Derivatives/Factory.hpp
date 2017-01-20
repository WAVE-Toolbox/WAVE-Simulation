

#pragma once

#include "Derivatives.hpp"
#include "FDTD2D.hpp"
#include "FDTD3D.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        namespace Derivatives {
            
            
            template<typename ValueType>
            class Factory
            {
                
            public:
                
                typedef typename Derivatives<ValueType>::DerivativesPtr  DerivativesPtr;
                
                Factory() = delete;
                Factory(Factory const&) = delete;
                void operator=(Factory const&) = delete;
                
                static DerivativesPtr Create( std::string const& dimension ) {
                    
                    // Assert correctness of input values
                    SCAI_ASSERT_ERROR( dimension.compare("3D")==0 || dimension.compare("2D")==0 , "Unkown dimension" );
                    
                    if( dimension.compare("2D") == 0 ){
                        return DerivativesPtr(new FDTD2D<ValueType>);
                    }
                    if( dimension.compare("3D") == 0 ){
                        return DerivativesPtr(new FDTD3D<ValueType>);
                    }
                    
                    COMMON_THROWEXCEPTION("Reached end of factory without match");
                    return nullptr;
                }
                
                
            };
        }
    }
}
