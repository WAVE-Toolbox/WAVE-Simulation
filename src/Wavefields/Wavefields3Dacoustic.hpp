/*! \brief The class Wavefields3Dacoustic holds the wavefields for 3D acoustic simulation
 *
 * Wavefields implements some methods, which are requiered by all derived classes.
 * As this class is an abstract class, all methods are protected.
 */


#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

template<typename ValueType>
class Wavefields3Dacoustic : public Wavefields<ValueType>
{
    
public:
    
    //! Default constructor
    Wavefields3Dacoustic(){};
    
    //! Default destructor
    ~Wavefields3Dacoustic(){};
    
    Wavefields3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    
    void reset();
    
    /* Getter routines for required wavefields */
    lama::DenseVector<ValueType>* getVX();
    lama::DenseVector<ValueType>* getVY();
    lama::DenseVector<ValueType>* getVZ();
    lama::DenseVector<ValueType>* getP();
    
    /* Getter routines for non-required wavefields: Will throw an error */
    lama::DenseVector<ValueType>* getSxx();
    lama::DenseVector<ValueType>* getSyy();
    lama::DenseVector<ValueType>* getSzz();
    lama::DenseVector<ValueType>* getSyz();
    lama::DenseVector<ValueType>* getSxz();
    lama::DenseVector<ValueType>* getSxy();

private:
    
    lama::DenseVector<ValueType> VX; //!< Wavefield for velocity in x
    lama::DenseVector<ValueType> VY; //!< Wavefield for velocity in y
    lama::DenseVector<ValueType> VZ; //!< Wavefield for velocity in z
    lama::DenseVector<ValueType> P; //!< Wavefield for pressure
    
};




/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D acoustic wavefields
 *
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
Wavefields3Dacoustic<ValueType>::Wavefields3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(VX,ctx,dist);
    this->initWavefield(VY,ctx,dist);
    this->initWavefield(VZ,ctx,dist);
    this->initWavefield(P,ctx,dist);
}


/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void Wavefields3Dacoustic<ValueType>::reset()
{
    this->resetWavefield(VX);
    this->resetWavefield(VY);
    this->resetWavefield(VZ);
    this->resetWavefield(P);
}


//! \brief Getter routine for vX wavefield
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getVX(){
    return(&VX);
}

//! \brief Getter routine for vY wavefield
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getVY(){
    return(&VY);
}

//! \brief Getter routine for vZ wavefield
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getVZ(){
    return(&VZ);
}

//! \brief Getter routine for p wavefield
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getP(){
    return(&P);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSxx(){
    COMMON_THROWEXCEPTION("There is no Sxx wavefield in the 3D acoustic case.")
    return(NULL);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSyy(){
    COMMON_THROWEXCEPTION("There is no Syy wavefield in the 3D acoustic case.")
    return(NULL);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSzz(){
    COMMON_THROWEXCEPTION("There is no Szz wavefield in the 3D acoustic case.")
    return(NULL);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSyz(){
    COMMON_THROWEXCEPTION("There is no Syz wavefield in the 3D acoustic case.")
    return(NULL);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSxz(){
    COMMON_THROWEXCEPTION("There is no Sxz wavefield in the 3D acoustic case.")
    return(NULL);
}

//! \brief Not valid in the 3D acoustic case
template<typename ValueType>
lama::DenseVector<ValueType>* Wavefields3Dacoustic<ValueType>::getSxy(){
    COMMON_THROWEXCEPTION("There is no Syx wavefield in the 3D acoustic case.")
    return(NULL);
}
