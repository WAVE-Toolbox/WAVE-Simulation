

#pragma once

#include "Modelparameter.hpp"
#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "Viscoelastic.hpp"

namespace KITGPI {
    
    namespace Modelparameter{
        
        
        template<typename ValueType>
        class Factory
        {
            
        public:
            
            typedef typename Modelparameter<ValueType>::ModelparameterPtr  ModelparameterPtr;
            
            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            static ModelparameterPtr Create( std::string const& type ) {
                
                // Assert correctness of input values
                SCAI_ASSERT_ERROR(type.compare("acoustic")==0 || type.compare("elastic")==0 || type.compare("visco")==0, "Unkown type" );
                
                if( type.compare("acoustic") == 0 ){
                    return ModelparameterPtr(new Acoustic<ValueType>);
                }
                if( type.compare("elastic") == 0 ){
                    return ModelparameterPtr(new Elastic<ValueType>);
                }
                if( type.compare("visco") == 0 ){
                    return ModelparameterPtr(new Viscoelastic<ValueType>);
                }
                
                COMMON_THROWEXCEPTION("Reached end of factory without match");
                return nullptr;
            }
            
            
        };
    }
}
