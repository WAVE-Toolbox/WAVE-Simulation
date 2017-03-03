

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

        /*! \brief The class FD3Delastic holds the wavefields for 3D elastic simulation
         *
         */
        template <typename ValueType>
        class FD3Delastic : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD3Delastic(){};

            //! Default destructor
            ~FD3Delastic(){};

            explicit FD3Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void reset() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getP() override;
            scai::lama::DenseVector<ValueType> &getRxx() override;
            scai::lama::DenseVector<ValueType> &getRyy() override;
            scai::lama::DenseVector<ValueType> &getRzz() override;
            scai::lama::DenseVector<ValueType> &getRyz() override;
            scai::lama::DenseVector<ValueType> &getRxz() override;
            scai::lama::DenseVector<ValueType> &getRxy() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;
            
            /* Overloading Operators */
            KITGPI::Wavefields::FD3Delastic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD3Delastic<ValueType> operator*=(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD3Delastic<ValueType> operator*(KITGPI::Wavefields::FD3Delastic<ValueType> rhs);
            KITGPI::Wavefields::FD3Delastic<ValueType> operator*=(KITGPI::Wavefields::FD3Delastic<ValueType> rhs);

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

            /* non-required wavefields */
            using Wavefields<ValueType>::P;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;
        };
    }
}
