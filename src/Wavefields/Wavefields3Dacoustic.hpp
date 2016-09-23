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
class Wavefields3Dacoustic : private Wavefields<ValueType>
{
    
public:
    
    void reset();
    
    Wavefields3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    //! Default deconstructor
    ~Wavefields3Dacoustic(){};
    
    lama::DenseVector<ValueType> vX; //!< Wavefield for velocity in x
    lama::DenseVector<ValueType> vY; //!< Wavefield for velocity in y
    lama::DenseVector<ValueType> vZ; //!< Wavefield for velocity in z
    lama::DenseVector<ValueType> p; //!< Wavefield for pressure
    
};

/*! \brief Constructor which will set context, allocate and set the wavefields to zero.
 *
 * Initialisation of 3D acoustic wavefields
 /param ctx Context
 /param dist Distribution
 */
template<typename ValueType>
Wavefields3Dacoustic<ValueType>::Wavefields3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    this->initWavefield(vX,ctx,dist);
    this->initWavefield(vY,ctx,dist);
    this->initWavefield(vZ,ctx,dist);
    this->initWavefield(p,ctx,dist);
}


/*! \brief Set all wavefields to zero.
 */
template<typename ValueType>
void Wavefields3Dacoustic<ValueType>::reset()
{
    this->resetWavefield(vX);
    this->resetWavefield(vY);
    this->resetWavefield(vZ);
    this->resetWavefield(p);
}
