

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
            
            FD3Delastic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset();
            
            /* Getter routines for required wavefields */
            lama::DenseVector<ValueType>& getVX();
            lama::DenseVector<ValueType>& getVY();
            lama::DenseVector<ValueType>& getVZ();
            lama::DenseVector<ValueType>& getSxx();
            lama::DenseVector<ValueType>& getSyy();
            lama::DenseVector<ValueType>& getSzz();
            lama::DenseVector<ValueType>& getSyz();
            lama::DenseVector<ValueType>& getSxz();
            lama::DenseVector<ValueType>& getSxy();
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getP();
            
            
        private:
            
            /* required wavefields */
            lama::DenseVector<ValueType> VX; //!< Wavefield for velocity in x
            lama::DenseVector<ValueType> VY; //!< Wavefield for velocity in y
            lama::DenseVector<ValueType> VZ; //!< Wavefield for velocity in z
            lama::DenseVector<ValueType> Sxx; //!< Wavefield (not-used here)
            lama::DenseVector<ValueType> Syy; //!< Wavefield (not-used here)
            lama::DenseVector<ValueType> Szz; //!< Wavefield (not-used here)
            lama::DenseVector<ValueType> Syz; //!< Wavefield (not-used here)
            lama::DenseVector<ValueType> Sxz; //!< Wavefield (not-used here)
            lama::DenseVector<ValueType> Sxy; //!< Wavefield (not-used here)
            
            /* non-required wavefields */
            lama::DenseVector<ValueType> P; //!< Wavefield (not-used here)

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


//! \brief Getter routine for vX wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getVX(){
    return(VX);
}

//! \brief Getter routine for vY wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getVY(){
    return(VY);
}

//! \brief Getter routine for vZ wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getVZ(){
    return(VZ);
}

//! \brief Not valid in the 3D elastic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getP(){
    COMMON_THROWEXCEPTION("There is no p wavefield in the 3D elastic case.")
    return(P);
}

//! \brief Getter routine for Sxx wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSxx(){
    return(Sxx);
}

//! \brief Getter routine for Syy wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSyy(){
    return(Syy);
}

//! \brief Getter routine for Szz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSzz(){
    return(Szz);
}

//! \brief Getter routine for Syz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSyz(){
    return(Syz);
}

//! \brief Getter routine for Sxz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSxz(){
    return(Sxz);
}

//! \brief Getter routine for Sxy wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Delastic<ValueType>::getSxy(){
    return(Sxy);
}
