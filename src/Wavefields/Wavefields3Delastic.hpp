

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD3Delastic holds the wavefields for 3D elastic simulation
         *
         */
        template<typename ValueType>
        class FD3Delastic : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD3Delastic(){};
            
            //! Default destructor
            ~FD3Delastic(){};
            
            explicit FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset();
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getP();
            lama::DenseVector<ValueType>& getRxx();
            lama::DenseVector<ValueType>& getRyy();
            lama::DenseVector<ValueType>& getRzz();
            lama::DenseVector<ValueType>& getRyz();
            lama::DenseVector<ValueType>& getRxz();
            lama::DenseVector<ValueType>& getRxy();
            
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
            
            /* non-required wavefields */
            using Wavefields<ValueType>::P; //!< Wavefield
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;
        };
    }
}


/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D elastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD3Delastic<ValueType>::FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
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
    
}

/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD3Delastic<ValueType>::reset()
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
}


//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getP(){
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D elastic case.")
    return(P);
}


//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRxx(){
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 3D elastic case.")
    return(Rxx);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRyy(){
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 3D elastic case.")
    return(Ryy);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRzz(){
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 3D elastic case.")
    return(Rzz);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRyz(){
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 3D elastic case.")
    return(Ryz);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRxz(){
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 3D elastic case.")
    return(Rxz);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getRxy(){
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 3D elastic case.")
    return(Rxy);
}
