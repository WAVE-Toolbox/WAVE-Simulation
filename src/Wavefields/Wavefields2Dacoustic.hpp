

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI {
    
    namespace Wavefields {
        
        /*! \brief The class FD2Dacoustic holds the wavefields for 2D acoustic simulation
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template<typename ValueType>
        class FD2Dacoustic : public Wavefields<ValueType>
        {
            
        public:
            
            //! Default constructor
            FD2Dacoustic(){};
            
            //! Default destructor
            ~FD2Dacoustic(){};
            
            FD2Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void reset();
            
            /* Getter routines for non-required wavefields: Will throw an error */
            lama::DenseVector<ValueType>& getVZ();
            lama::DenseVector<ValueType>& getSxx();
            lama::DenseVector<ValueType>& getSyy();
            lama::DenseVector<ValueType>& getSzz();
            lama::DenseVector<ValueType>& getSyz();
            lama::DenseVector<ValueType>& getSxz();
            lama::DenseVector<ValueType>& getSxy();
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
            using Wavefields<ValueType>::P;
            
            /* non-required wavefields */
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
        };
    }
}


/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 2D acoustic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
KITGPI::Wavefields::FD2Dacoustic<ValueType>::FD2Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(P,ctx,dist);
}


/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::FD2Dacoustic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(P);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getVZ(){
    COMMON_THROWEXCEPTION("There is no Vz wavefield in the 2D acoustic case.")
    return(VZ);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxx(){
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 2D acoustic case.")
    return(Sxx);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSyy(){
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 2D acoustic case.")
    return(Syy);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSzz(){
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 2D acoustic case.")
    return(Szz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSyz(){
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 2D acoustic case.")
    return(Syz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxz(){
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 2D acoustic case.")
    return(Sxz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getSxy(){
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 2D acoustic case.")
    return(Sxy);
}


//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxx(){
    COMMON_THROWEXCEPTION("There is no Rxx wavefield in the 2D acoustic case.")
    return(Rxx);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRyy(){
    COMMON_THROWEXCEPTION("There is no Ryy wavefield in the 2D acoustic case.")
    return(Ryy);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRzz(){
    COMMON_THROWEXCEPTION("There is no Rzz wavefield in the 2D acoustic case.")
    return(Rzz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRyz(){
    COMMON_THROWEXCEPTION("There is no Ryz wavefield in the 2D acoustic case.")
    return(Ryz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxz(){
    COMMON_THROWEXCEPTION("There is no Rxz wavefield in the 2D acoustic case.")
    return(Rxz);
}

//! \brief Not valid in the 2D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::FD2Dacoustic<ValueType>::getRxy(){
    COMMON_THROWEXCEPTION("There is no Rxy wavefield in the 2D acoustic case.")
    return(Rxy);
}
