

#pragma once

#include "TMEM.hpp"
#include "EMEM.hpp"
#include "ViscoTMEM.hpp"
#include "ViscoEMEM.hpp"
#include "Modelparameter.hpp"
#include <string>

namespace KITGPI
{

    namespace Modelparameter
    {

        //! \brief Factory class.
        template <typename ValueType>
        class FactoryEM
        {

          public:
            //! \brief Declare Modelparameter pointer
            typedef typename ModelparameterEM<ValueType>::ModelparameterPtr ModelparameterPtr;

            FactoryEM() = delete;
            FactoryEM(FactoryEM const &) = delete;
            void operator=(FactoryEM const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param type Simmulation type (tmem, emem, viscoemem)
             */
            static ModelparameterPtr Create(std::string type);
        };
    }
}
