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

        /*! \brief The class FD2Dviscoelastic holds the wavefields for 2D viscoelastic simulation
         *
         */
        template <typename ValueType>
        class FD2Dviscoelastic : public WavefieldsSeismic<ValueType>
        {

          public:
            //! Default constructor
            FD2Dviscoelastic():EquationType("viscoelastic"),NumDimension(2)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD2Dviscoelastic(){};

            explicit FD2Dviscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefP() override;
            scai::lama::DenseVector<ValueType> &getRefVZ() override;
            scai::lama::DenseVector<ValueType> &getRefSzz() override;
            scai::lama::DenseVector<ValueType> &getRefSyz() override;
            scai::lama::DenseVector<ValueType> &getRefSxz() override;
            scai::lama::DenseVector<ValueType> &getRefRzz() override;
            scai::lama::DenseVector<ValueType> &getRefRyz() override;
            scai::lama::DenseVector<ValueType> &getRefRxz() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dviscoelastic<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoelastic<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoelastic<ValueType> operator*(KITGPI::Wavefields::FD2Dviscoelastic<ValueType> rhs);
            KITGPI::Wavefields::FD2Dviscoelastic<ValueType> operator*=(KITGPI::Wavefields::FD2Dviscoelastic<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &SWaveModulus);

            std::string EquationType;
            int NumDimension;
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;

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

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
