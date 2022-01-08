#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "WavefieldsSeismic.hpp"

namespace KITGPI
{

    namespace Wavefields
    {

        /*! \brief The class FD3Dacoustic holds the wavefields for 3D acoustic simulation
         *
         */
        template <typename ValueType>
        class FD3Dacoustic : public WavefieldsSeismic<ValueType>
        {

          public:
            //! Default constructor
            FD3Dacoustic():EquationType("acoustic"),NumDimension(3)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD3Dacoustic(){};

            explicit FD3Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefSxx() override;
            scai::lama::DenseVector<ValueType> &getRefSyy() override;
            scai::lama::DenseVector<ValueType> &getRefSzz() override;
            scai::lama::DenseVector<ValueType> &getRefSyz() override;
            scai::lama::DenseVector<ValueType> &getRefSxz() override;
            scai::lama::DenseVector<ValueType> &getRefSxy() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRxx() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRyy() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRzz() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRyz() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRxz() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRxy() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            ValueType estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*(KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs);
            KITGPI::Wavefields::FD3Dacoustic<ValueType> operator*=(KITGPI::Wavefields::FD3Dacoustic<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposeWavefieldType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
              std::string EquationType;
              int NumDimension;
            using Wavefields<ValueType>::equationType; 
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::numRelaxationMechanisms;

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

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
