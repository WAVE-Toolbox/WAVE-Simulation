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

        /*! \brief The class FD3Dviscoelastic holds the wavefields for 3D viscoelastic simulation
         *
         */
        template <typename ValueType>
        class FD3Dviscoelastic : public WavefieldsSeismic<ValueType>
        {

          public:
            //! Default constructor
            FD3Dviscoelastic():EquationType("viscoelastic"),NumDimension(3)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD3Dviscoelastic(){};

            explicit FD3Dviscoelastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefP() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            ValueType estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Dviscoelastic<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD3Dviscoelastic<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD3Dviscoelastic<ValueType> operator*(KITGPI::Wavefields::FD3Dviscoelastic<ValueType> rhs);
            KITGPI::Wavefields::FD3Dviscoelastic<ValueType> operator*=(KITGPI::Wavefields::FD3Dviscoelastic<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void timesAssign(scai::lama::DenseVector<ValueType> rhs);
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &SWaveModulus);

            std::string EquationType;
            int NumDimension;
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;
            using Wavefields<ValueType>::numRelaxationMechanisms;

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

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
