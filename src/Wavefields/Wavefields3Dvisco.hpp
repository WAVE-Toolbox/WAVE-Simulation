

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"

namespace KITGPI
{

    namespace Wavefields
    {

        /*! \brief The class FD3Dvisco holds the wavefields for 3D visco elastic simulation
         *
         */
        template <typename ValueType>
        class FD3Dvisco : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD3Dvisco(){};

            //! Default destructor
            ~FD3Dvisco(){};

            explicit FD3Dvisco(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefP() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Dvisco<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD3Dvisco<ValueType> operator*=(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD3Dvisco<ValueType> operator*(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs);
            KITGPI::Wavefields::FD3Dvisco<ValueType> operator*=(KITGPI::Wavefields::FD3Dvisco<ValueType> rhs);

            void write(std::string type, IndexType t, IndexType partitionedOut) override;
            void writeSnapshot(IndexType t, IndexType partitionedOut);

	    void minusAssign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
	    void timesAssign(ValueType rhs);
	    
          private:
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
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

            /* non-required wavefields */
            using Wavefields<ValueType>::P; //!< Wavefield

            std::string type = "Visco3D";
        };
    }
}
