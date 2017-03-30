

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

        /*! \brief The class FD3Dacoustic holds the wavefields for 3D acoustic simulation
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class FD3Dacoustic : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD3Dacoustic(){};

            //! Default destructor
            ~FD3Dacoustic(){};

            explicit FD3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void reset() override;

            /* Getter routines for non-required wavefields: Will throw an error */
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

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);

          private:
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::P;

            /* non-required wavefields */
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

            std::string type = "Acoustic3D";
        };
    }
}
