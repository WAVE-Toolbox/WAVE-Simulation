

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

        /*! \brief The class FD2Dsh holds the wavefields for 2D sh simulation
         *
         */
        template <typename ValueType>
        class FD2Dsh : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD2Dsh(){};

            //! Default destructor
            ~FD2Dsh(){};

            explicit FD2Dsh(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void reset() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getP() override;
            scai::lama::DenseVector<ValueType> &getRxx() override;
            scai::lama::DenseVector<ValueType> &getRyy() override;
            scai::lama::DenseVector<ValueType> &getRzz() override;
            scai::lama::DenseVector<ValueType> &getRyz() override;
            scai::lama::DenseVector<ValueType> &getRxz() override;
            scai::lama::DenseVector<ValueType> &getRxy() override;
            scai::lama::DenseVector<ValueType> &getVX() override;
            scai::lama::DenseVector<ValueType> &getVY() override;
            scai::lama::DenseVector<ValueType> &getSxx() override;
            scai::lama::DenseVector<ValueType> &getSyy() override;
            scai::lama::DenseVector<ValueType> &getSzz() override;
            scai::lama::DenseVector<ValueType> &getSxy() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dsh<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD2Dsh<ValueType> operator*=(scai::lama::Scalar rhs);
            KITGPI::Wavefields::FD2Dsh<ValueType> operator*(KITGPI::Wavefields::FD2Dsh<ValueType> rhs);
            KITGPI::Wavefields::FD2Dsh<ValueType> operator*=(KITGPI::Wavefields::FD2Dsh<ValueType> rhs);

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);

          private:
            /* required wavefields */
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::Syz;

            /* non-required wavefields */
            using Wavefields<ValueType>::P;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Sxy;
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;

            std::string type = "SH2D";
        };
    }
}
