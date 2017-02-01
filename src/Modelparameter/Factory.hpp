

#pragma once

#include <string>
#include "Modelparameter.hpp"
#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "Viscoelastic.hpp"

namespace KITGPI {
    
    namespace Modelparameter{
        
        //! \brief Factory class.
        template<typename ValueType>
        class Factory
        {
            
        public:
            
            //! \brief Declare Modelparameter pointer
            typedef typename Modelparameter<ValueType>::ModelparameterPtr  ModelparameterPtr;
            
            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            /*! \brief Create the right simmulation with factory methode.
             *
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ModelparameterPtr Create( std::string type ) {
                
                // transform to lower cases
                std::transform(type.begin(), type.end(), type.begin(), ::tolower);
                
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
