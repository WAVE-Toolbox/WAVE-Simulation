#pragma once

#include "Acoustic.hpp"
#include "Elastic.hpp"
#include "SH.hpp"
#include "Viscoelastic.hpp"
#include "Modelparameter.hpp"
#include "ModelparameterSeismic.hpp"
#include "../ModelparameterEM/TMEM.hpp"
#include "../ModelparameterEM/EMEM.hpp"
#include "../ModelparameterEM/ViscoTMEM.hpp"
#include "../ModelparameterEM/ViscoEMEM.hpp"
#include "../ModelparameterEM/ModelparameterEM.hpp"
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
             \param type Simmulation type (acoustic, elsstic, viscoelastic, sh)
             */
            static ModelparameterPtr Create(std::string type);
        };
    }
}
