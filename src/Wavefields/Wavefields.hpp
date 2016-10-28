

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>


namespace KITGPI {
    
    //! \brief Wavefields namespace
    namespace Wavefields {
        
        /*! \brief Abstract class to handle the wavefields for the forward modelling.
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template<typename ValueType>
        class Wavefields
        {
            
        public:
            
            //! Default deconstructor
            ~Wavefields(){};
            
            //! Reset wavefields
            virtual void reset()=0;
            
            //! \brief Getter routine for vX wavefield
            virtual lama::DenseVector<ValueType>& getVX();
            //! \brief Getter routine for vY wavefield
            virtual lama::DenseVector<ValueType>& getVY();
            //! \brief Getter routine for vZ wavefield
            virtual lama::DenseVector<ValueType>& getVZ();
            
            //! \brief Getter routine for p wavefield
            virtual lama::DenseVector<ValueType>& getP();
            
            //! \brief Getter routine for sxx wavefield
            virtual lama::DenseVector<ValueType>& getSxx();
            //! \brief Getter routine for syy wavefield
            virtual lama::DenseVector<ValueType>& getSyy();
            //! \brief Getter routine for szz wavefield
            virtual lama::DenseVector<ValueType>& getSzz();
            
            //! \brief Getter routine for syx wavefield
            virtual lama::DenseVector<ValueType>& getSyz();
            //! \brief Getter routine for sxz wavefield
            virtual lama::DenseVector<ValueType>& getSxz();
            //! \brief Getter routine for sxy wavefield
            virtual lama::DenseVector<ValueType>& getSxy();
            
        protected:
            
            void resetWavefield(lama::DenseVector<ValueType>& vector);
            void initWavefield(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            lama::DenseVector<ValueType> VX; //!< Wavefield for velocity in x
            lama::DenseVector<ValueType> VY; //!< Wavefield for velocity in y
            lama::DenseVector<ValueType> VZ; //!< Wavefield for velocity in z
            lama::DenseVector<ValueType> Sxx; //!< Wavefield
            lama::DenseVector<ValueType> Syy; //!< Wavefield
            lama::DenseVector<ValueType> Szz; //!< Wavefield
            lama::DenseVector<ValueType> Syz; //!< Wavefield
            lama::DenseVector<ValueType> Sxz; //!< Wavefield
            lama::DenseVector<ValueType> Sxy; //!< Wavefield
            lama::DenseVector<ValueType> P; //!< Wavefield
            
        };
    }
}

/*! \brief Reset a single wavefield to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::resetWavefield(lama::DenseVector<ValueType>& vector)
{
    vector.assign(0.0);
}

/*! \brief Intitialisation of a single wavefield vector.
 *
 * This method will set the context, allocate the the wavefield and set the field to zero.
 */
template<typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::initWavefield(lama::DenseVector<ValueType>& vector,hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
    
    resetWavefield(vector);
}

//! \brief Getter routine for vX wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getVX(){
    return(VX);
}

//! \brief Getter routine for vY wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getVY(){
    return(VY);
}

//! \brief Getter routine for vZ wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getVZ(){
    return(VZ);
}


//! \brief Getter routine for Sxx wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSxx(){
    return(Sxx);
}

//! \brief Getter routine for Syy wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSyy(){
    return(Syy);
}

//! \brief Getter routine for Szz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSzz(){
    return(Szz);
}

//! \brief Getter routine for Syz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSyz(){
    return(Syz);
}

//! \brief Getter routine for Sxz wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSxz(){
    return(Sxz);
}

//! \brief Getter routine for Sxy wavefield
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getSxy(){
    return(Sxy);
}

//! \brief Getter routine for P
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getP(){
    return(P);
}
