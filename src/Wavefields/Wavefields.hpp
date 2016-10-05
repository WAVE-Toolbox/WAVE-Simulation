/*! \brief The class Wavefields is a abstract type that represents the wavefield for the forward modelling.
 *
 * Wavefields implements some methods, which are requiered by all derived classes.
 * As this class is an abstract class, all methods are protected.
 */

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

template<typename ValueType>
class Wavefields
{
    
protected:
    
    //! Default constructor
    Wavefields(){};
    //! Default deconstructor
    ~Wavefields(){};
    
    //! Reset wavefields
    virtual void reset()=0;
    
    void resetWavefield(lama::DenseVector<ValueType>& vector);
    void initWavefield(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    
};


/*! \brief Reset a single wavefield to zero.
 */
template<typename ValueType>
void Wavefields<ValueType>::resetWavefield(lama::DenseVector<ValueType>& vector)
{
    vector.assign(0.0);
}

/*! \brief Intitialisation of a single wavefield vector.
 *
 * This method will set the context, allocate the the wavefield and set the field to zero.
 */
template<typename ValueType>
void Wavefields<ValueType>::initWavefield(lama::DenseVector<ValueType>& vector,hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
    
    resetWavefield(vector);
}

