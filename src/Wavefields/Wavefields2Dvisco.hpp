

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

        /*! \brief The class FD2Dvisco holds the wavefields for 2D visco elastic simulation
         *
         */
        template <typename ValueType>
        class FD2Dvisco : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD2Dvisco(){};

            //! Default destructor
            ~FD2Dvisco(){};

            explicit FD2Dvisco(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void reset() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getP() override;
            scai::lama::DenseVector<ValueType> &getVZ() override;
            scai::lama::DenseVector<ValueType> &getSzz() override;
            scai::lama::DenseVector<ValueType> &getSyz() override;
            scai::lama::DenseVector<ValueType> &getSxz() override;
            scai::lama::DenseVector<ValueType> &getRzz() override;
            scai::lama::DenseVector<ValueType> &getRyz() override;
            scai::lama::DenseVector<ValueType> &getRxz() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string type, IndexType t) override;
            void writeSnapshot(IndexType t);

          private:
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Sxy;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rxy;

            /* non-required wavefields */
            using Wavefields<ValueType>::P; //!< Wavefield
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rzz;

            std::string type = "Visco2D";
        };
    }
}
