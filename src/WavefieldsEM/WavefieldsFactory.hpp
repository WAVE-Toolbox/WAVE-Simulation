#pragma once

#include "Wavefields.hpp"
#include "Wavefields2Dtmem.hpp"
#include "Wavefields2Demem.hpp"
#include "Wavefields2Dviscotmem.hpp"
#include "Wavefields2Dviscoemem.hpp"
#include "Wavefields3Demem.hpp"
#include "Wavefields3Dviscoemem.hpp"

namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        //! \brief FactoryEM class.
        template <typename ValueType>
        class FactoryEM
        {

          public:
            //! \brief Declare Wavefield pointer
            typedef typename WavefieldsEM<ValueType>::WavefieldPtr WavefieldPtrEM;

            FactoryEM() = delete;
            FactoryEM(FactoryEM const &) = delete;
            void operator=(FactoryEM const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (tmem, emem, viscoemem)
             */
            static WavefieldPtrEM Create(std::string dimension, std::string type);
        };
    }
}
