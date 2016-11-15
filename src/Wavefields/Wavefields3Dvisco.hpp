

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD3Dvisco holds the wavefields for 3D visco elastic simulation
         *
         */
        template<typename ValueType>
        class FD3Dvisco : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD3Dvisco(){};
            
            //! Default destructor
            ~FD3Dvisco(){};
            
            FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset();
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getP();
            
        private:
            
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::Sxy;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;
            
            /* non-required wavefields */
            using Wavefields<ValueType>::P; //!< Wavefield

        };
    }
}


/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D viscoelastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD3Dvisco<ValueType>::FD3Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(VZ,ctx,dist);
    this->initWavefield(Sxx,ctx,dist);
    this->initWavefield(Syy,ctx,dist);
    this->initWavefield(Szz,ctx,dist);
    this->initWavefield(Syz,ctx,dist);
    this->initWavefield(Sxz,ctx,dist);
    this->initWavefield(Sxy,ctx,dist);
    this->initWavefield(Rxx,ctx,dist);
    this->initWavefield(Ryy,ctx,dist);
    this->initWavefield(Rzz,ctx,dist);
    this->initWavefield(Ryz,ctx,dist);
    this->initWavefield(Rxz,ctx,dist);
    this->initWavefield(Rxy,ctx,dist);

}

/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD3Dvisco<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Szz);
    this->resetWavefield(Syz);
    this->resetWavefield(Sxz);
    this->resetWavefield(Sxy);
    this->resetWavefield(Rxx);
    this->resetWavefield(Ryy);
    this->resetWavefield(Rzz);
    this->resetWavefield(Ryz);
    this->resetWavefield(Rxz);
    this->resetWavefield(Rxy);
}


//! \brief Not valid in the 3D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dvisco<ValueType>::getP(){
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D visco-elastic case.")
    return(P);
}

