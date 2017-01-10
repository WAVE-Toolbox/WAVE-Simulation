

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD2Dvisco holds the wavefields for 2D visco elastic simulation
         *
         */
        template<typename ValueType>
        class FD2Dvisco : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD2Dvisco(){};
            
            //! Default destructor
            ~FD2Dvisco(){};
            
            explicit FD2Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset() override;
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getP() override;
            lama::DenseVector<ValueType>& getVZ() override;
            lama::DenseVector<ValueType>& getSzz() override;
            lama::DenseVector<ValueType>& getSyz() override;
            lama::DenseVector<ValueType>& getSxz() override;
            lama::DenseVector<ValueType>& getRzz() override;
            lama::DenseVector<ValueType>& getRyz() override;
            lama::DenseVector<ValueType>& getRxz() override;
            
            hmemo::ContextPtr getContextPtr() override;
            
        private:
            
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Sxy;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rxy;
            
            /* non-required wavefields */
            using Wavefields<ValueType>::P; //!< Wavefield
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rzz;


        };
    }
}


/*! \brief Returns hmemo::ContextPtr from this wavefields
 */
template<typename ValueType>
hmemo::ContextPtr KITGPI::Wavefields::FD2Dvisco<ValueType>::getContextPtr()
{
    return(VX.getContextPtr());
}

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D viscoelastic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD2Dvisco<ValueType>::FD2Dvisco(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(Sxx,ctx,dist);
    this->initWavefield(Syy,ctx,dist);
    this->initWavefield(Sxy,ctx,dist);
    this->initWavefield(Rxx,ctx,dist);
    this->initWavefield(Ryy,ctx,dist);
    this->initWavefield(Rxy,ctx,dist);

}

/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD2Dvisco<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(Sxx);
    this->resetWavefield(Syy);
    this->resetWavefield(Sxy);
    this->resetWavefield(Rxx);
    this->resetWavefield(Ryy);
    this->resetWavefield(Rxy);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getRzz(){
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D visco-elastic case.")
    return(Rzz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getRyz(){
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D visco-elastic case.")
    return(Ryz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getRxz(){
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D visco-elastic case.")
    return(Rxz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getSzz(){
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D visco-elastic case.")
    return(Szz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getSyz(){
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D visco-elastic case.")
    return(Syz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getSxz(){
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D visco-elastic case.")
    return(Sxz);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getVZ(){
    COMMON_THROWEXCEPTION("There is no VZ wavefield in the 2D visco-elastic case.")
    return(VZ);
}

//! \brief Not valid in the 2D visco-elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dvisco<ValueType>::getP(){
    COMMON_THROWEXCEPTION("There is no p wavefield in the 2D visco-elastic case.")
    return(P);
}

