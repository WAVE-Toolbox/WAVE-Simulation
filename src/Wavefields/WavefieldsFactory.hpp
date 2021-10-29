#pragma once

#include "Wavefields.hpp"
#include "WavefieldsSeismic.hpp"
#include "Wavefields2Dacoustic.hpp"
#include "Wavefields2Delastic.hpp"
#include "Wavefields2Dsh.hpp"
#include "Wavefields2Dviscosh.hpp"
#include "Wavefields2Dviscoelastic.hpp"
#include "Wavefields3Dacoustic.hpp"
#include "Wavefields3Delastic.hpp"
#include "Wavefields3Dviscoelastic.hpp"
#include "../WavefieldsEM/WavefieldsEM.hpp"
#include "../WavefieldsEM/Wavefields2Dtmem.hpp"
#include "../WavefieldsEM/Wavefields2Demem.hpp"
#include "../WavefieldsEM/Wavefields2Dviscotmem.hpp"
#include "../WavefieldsEM/Wavefields2Dviscoemem.hpp"
#include "../WavefieldsEM/Wavefields3Demem.hpp"
#include "../WavefieldsEM/Wavefields3Dviscoemem.hpp"

namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        //! \brief Factory class.
        template <typename ValueType>
        class Factory
        {

          public:
            //! \brief Declare Wavefield pointer
            typedef typename Wavefields<ValueType>::WavefieldPtr WavefieldPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (acoustic, elsstic, viscoelastic, sh)
             */
            static WavefieldPtr Create(std::string dimension, std::string type);
        };
    }
}
