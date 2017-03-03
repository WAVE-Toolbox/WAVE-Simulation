

#pragma once

#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "Modelparameter.hpp"
#include "Viscoelastic.hpp"
#include <string>

namespace KITGPI
{

    namespace Modelparameter
    {

        //! \brief Factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare Modelparameter pointer
            typedef typename Modelparameter<ValueType>::ModelparameterPtr ModelparameterPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ModelparameterPtr Create(std::string type);
        };
    }
}
