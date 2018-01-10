

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        /*! \brief Abstract class to handle the wavefields for the forward modelling.
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class Wavefields
        {

          public:
            //! Default constructor
            Wavefields(){};
            //! Default destructor
            ~Wavefields(){};

            //! \brief Declare Wavefield pointer
            typedef std::shared_ptr<Wavefields<ValueType>> WavefieldPtr;

            //! Reset wavefields
            virtual void reset() = 0;

            virtual scai::lama::DenseVector<ValueType> &getVX();
            virtual scai::lama::DenseVector<ValueType> &getVY();
            virtual scai::lama::DenseVector<ValueType> &getVZ();
            virtual scai::lama::DenseVector<ValueType> &getP();

            virtual scai::lama::DenseVector<ValueType> &getSxx();
            virtual scai::lama::DenseVector<ValueType> &getSyy();
            virtual scai::lama::DenseVector<ValueType> &getSzz();
            virtual scai::lama::DenseVector<ValueType> &getSyz();
            virtual scai::lama::DenseVector<ValueType> &getSxz();
            virtual scai::lama::DenseVector<ValueType> &getSxy();

            virtual scai::lama::DenseVector<ValueType> &getRxx();
            virtual scai::lama::DenseVector<ValueType> &getRyy();
            virtual scai::lama::DenseVector<ValueType> &getRzz();
            virtual scai::lama::DenseVector<ValueType> &getRyz();
            virtual scai::lama::DenseVector<ValueType> &getRxz();
            virtual scai::lama::DenseVector<ValueType> &getRxy();

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void write(std::string type, IndexType t) = 0;

            //! Operator overloading
            virtual void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;

            KITGPI::Wavefields::Wavefields<ValueType> &operator=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator-=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator+=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, IndexType t);

            scai::lama::DenseVector<ValueType> VX;  //!< Wavefield for velocity in x
            scai::lama::DenseVector<ValueType> VY;  //!< Wavefield for velocity in y
            scai::lama::DenseVector<ValueType> VZ;  //!< Wavefield for velocity in z
            scai::lama::DenseVector<ValueType> Sxx; //!< Wavefield
            scai::lama::DenseVector<ValueType> Syy; //!< Wavefield
            scai::lama::DenseVector<ValueType> Szz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Syz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Sxz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Sxy; //!< Wavefield
            scai::lama::DenseVector<ValueType> P;   //!< Wavefield

            scai::lama::DenseVector<ValueType> Rxx; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryy; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rzz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxy; //!< Relaxation parameter
        };
    }
}
