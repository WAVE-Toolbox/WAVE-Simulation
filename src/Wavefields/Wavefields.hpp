

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
            
            virtual lama::DenseVector<ValueType>& getVX();
            virtual lama::DenseVector<ValueType>& getVY();
            virtual lama::DenseVector<ValueType>& getVZ();
            virtual lama::DenseVector<ValueType>& getP();
            
            virtual lama::DenseVector<ValueType>& getSxx();
            virtual lama::DenseVector<ValueType>& getSyy();
            virtual lama::DenseVector<ValueType>& getSzz();
            virtual lama::DenseVector<ValueType>& getSyz();
            virtual lama::DenseVector<ValueType>& getSxz();
            virtual lama::DenseVector<ValueType>& getSxy();
            
            virtual lama::DenseVector<ValueType>& getRxx();
            virtual lama::DenseVector<ValueType>& getRyy();
            virtual lama::DenseVector<ValueType>& getRzz();
            virtual lama::DenseVector<ValueType>& getRyz();
            virtual lama::DenseVector<ValueType>& getRxz();
            virtual lama::DenseVector<ValueType>& getRxy();
            
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
            
            lama::DenseVector<ValueType> Rxx; //!< Relaxation parameter
            lama::DenseVector<ValueType> Ryy; //!< Relaxation parameter
            lama::DenseVector<ValueType> Rzz; //!< Relaxation parameter
            lama::DenseVector<ValueType> Ryz; //!< Relaxation parameter
            lama::DenseVector<ValueType> Rxz; //!< Relaxation parameter
            lama::DenseVector<ValueType> Rxy; //!< Relaxation parameter
            
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

//! \brief Getter routine for Rxx Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRxx(){
    return(Rxx);
}

//! \brief Getter routine for Ryy Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRyy(){
    return(Ryy);
}

//! \brief Getter routine for Rzz Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRzz(){
    return(Rzz);
}

//! \brief Getter routine for Ryz Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRyz(){
    return(Ryz);
}

//! \brief Getter routine for Rxz Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRxz(){
    return(Rxz);
}

//! \brief Getter routine for Rxy Relaxation parameter
template<typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Wavefields::Wavefields<ValueType>::getRxy(){
    return(Rxy);
}
