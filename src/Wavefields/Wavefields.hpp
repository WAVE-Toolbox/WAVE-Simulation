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
    
public:
    
    //! Default deconstructor
    ~Wavefields(){};
    
    //! Reset wavefields
    virtual void reset()=0;
    
    //! \brief Getter routine for vX wavefield
    virtual lama::DenseVector<ValueType>* getVX()=0;
    //! \brief Getter routine for vY wavefield
    virtual lama::DenseVector<ValueType>* getVY()=0;
    //! \brief Getter routine for vZ wavefield
    virtual lama::DenseVector<ValueType>* getVZ()=0;
    
    //! \brief Getter routine for p wavefield
    virtual lama::DenseVector<ValueType>* getP()=0;
    
    //! \brief Getter routine for sxx wavefield
    virtual lama::DenseVector<ValueType>* getSxx()=0;
    //! \brief Getter routine for syy wavefield
    virtual lama::DenseVector<ValueType>* getSyy()=0;
    //! \brief Getter routine for szz wavefield
    virtual lama::DenseVector<ValueType>* getSzz()=0;
    
    //! \brief Getter routine for syx wavefield
    virtual lama::DenseVector<ValueType>* getSyz()=0;
    //! \brief Getter routine for sxz wavefield
    virtual lama::DenseVector<ValueType>* getSxz()=0;
    //! \brief Getter routine for sxy wavefield
    virtual lama::DenseVector<ValueType>* getSxy()=0;
    
protected:
    
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

