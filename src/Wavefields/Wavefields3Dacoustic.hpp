

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD3Dacoustic holds the wavefields for 3D acoustic simulation
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template<typename ValueType>
        class FD3Dacoustic : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD3Dacoustic(){};
            
            //! Default destructor
            ~FD3Dacoustic(){};
            
            FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset();
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getSxx();
            lama::DenseVector<ValueType>& getSyy();
            lama::DenseVector<ValueType>& getSzz();
            lama::DenseVector<ValueType>& getSyz();
            lama::DenseVector<ValueType>& getSxz();
            lama::DenseVector<ValueType>& getSxy();
            
        private:
            
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::P;
            
            /* non-required wavefields */
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::Sxy;
        };
    }
}


/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D acoustic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD3Dacoustic<ValueType>::FD3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(VZ,ctx,dist);
    this->initWavefield(P,ctx,dist);
}


/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD3Dacoustic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(P);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxx(){
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 3D acoustic case.")
    return(Sxx);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSyy(){
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 3D acoustic case.")
    return(Syy);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSzz(){
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 3D acoustic case.")
    return(Szz);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSyz(){
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 3D acoustic case.")
    return(Syz);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxz(){
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 3D acoustic case.")
    return(Sxz);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD3Dacoustic<ValueType>::getSxy(){
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 3D acoustic case.")
    return(Sxy);
}
