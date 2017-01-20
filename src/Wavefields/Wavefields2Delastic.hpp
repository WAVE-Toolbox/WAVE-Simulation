

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD2Delastic holds the wavefields for 2D elastic simulation
         *
         */
        template<typename ValueType>
        class FD2Delastic : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD2Delastic(){};
            
            //! Default destructor
            ~FD2Delastic(){};
            
            explicit FD2Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset() override;
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getP() override;
            lama::DenseVector<ValueType>& getRxx() override;
            lama::DenseVector<ValueType>& getRyy() override;
            lama::DenseVector<ValueType>& getRzz() override;
            lama::DenseVector<ValueType>& getRyz() override;
            lama::DenseVector<ValueType>& getRxz() override;
            lama::DenseVector<ValueType>& getRxy() override;
            lama::DenseVector<ValueType>& getVZ() override;
            lama::DenseVector<ValueType>& getSzz() override;
            lama::DenseVector<ValueType>& getSyz() override;
            lama::DenseVector<ValueType>& getSxz() override;
            
            hmemo::ContextPtr getContextPtr() override;
            
            void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist) override;
            
        private:
            
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Sxy;
            
            /* non-required wavefields */
            using Wavefields<ValueType>::P;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;
        };
    }
}

/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template<typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Delastic<ValueType>::getContextPtr()
{
    return(VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D elastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD2Delastic<ValueType>::FD2Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    init(ctx,dist);
}

template<typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(Sxx,ctx,dist);
    this->initWavefield(Syy,ctx,dist);
    this->initWavefield(Sxy,ctx,dist);
}


/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD2Delastic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Sxy);
}


//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getSxz(){
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D elastic case.")
    return(Sxz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getSyz(){
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D elastic case.")
    return(Syz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getSzz(){
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D elastic case.")
    return(Szz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getVZ(){
    COMMON_THROWEXCEPTION("There is no VZ wavefield in the 2D elastic case.")
    return(VZ);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getP(){
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D elastic case.")
    return(P);
}


//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRxx(){
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D elastic case.")
    return(Rxx);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRyy(){
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D elastic case.")
    return(Ryy);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRzz(){
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D elastic case.")
    return(Rzz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRyz(){
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D elastic case.")
    return(Ryz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRxz(){
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D elastic case.")
    return(Rxz);
}

//! \brief Not valid in the 2D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Delastic<ValueType>::getRxy(){
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D elastic case.")
    return(Rxy);
}
