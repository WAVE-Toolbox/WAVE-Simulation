

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

        /*! \brief The class FD2Dacoustic holds the wavefields for 2D acoustic simulation
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class FD2Dacoustic : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD2Dacoustic(){};

            //! Default destructor
            ~FD2Dacoustic(){};

            explicit FD2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void reset() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getVZ() override;
            scai::lama::DenseVector<ValueType> &getSxx() override;
            scai::lama::DenseVector<ValueType> &getSyy() override;
            scai::lama::DenseVector<ValueType> &getSzz() override;
            scai::lama::DenseVector<ValueType> &getSyz() override;
            scai::lama::DenseVector<ValueType> &getSxz() override;
            scai::lama::DenseVector<ValueType> &getSxy() override;
            scai::lama::DenseVector<ValueType> &getRxx() override;
            scai::lama::DenseVector<ValueType> &getRyy() override;
            scai::lama::DenseVector<ValueType> &getRzz() override;
            scai::lama::DenseVector<ValueType> &getRyz() override;
            scai::lama::DenseVector<ValueType> &getRxz() override;
            scai::lama::DenseVector<ValueType> &getRxy() override;

            scai::hmemo::ContextPtr getContextPtr() override;
            
            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*=(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*=(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs);

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
