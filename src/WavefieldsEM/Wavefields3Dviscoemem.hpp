

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "Wavefields.hpp"

namespace KITGPI
{

    namespace Wavefields
    {

        /*! \brief The class FD3Dviscoemem holds the wavefields for 3D viscoemem emem simulation
         *
         */
        template <typename ValueType>
        class FD3Dviscoemem : public WavefieldsEM<ValueType>
        {

          public:
            //! Default constructor
            FD3Dviscoemem():EquationType("viscoemem"),NumDimension(3)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD3Dviscoemem(){};

            explicit FD3Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Dviscoemem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD3Dviscoemem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD3Dviscoemem<ValueType> operator*(KITGPI::Wavefields::FD3Dviscoemem<ValueType> rhs);
            KITGPI::Wavefields::FD3Dviscoemem<ValueType> operator*=(KITGPI::Wavefields::FD3Dviscoemem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void applyWavefieldTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &SWaveModulus);

            std::string EquationType;
            int NumDimension;
            using WavefieldsEM<ValueType>::numDimension;
            using WavefieldsEM<ValueType>::equationType;

            /* required wavefields */
            using WavefieldsEM<ValueType>::HX;
            using WavefieldsEM<ValueType>::HY;
            using WavefieldsEM<ValueType>::HZ;
            using WavefieldsEM<ValueType>::EX;
            using WavefieldsEM<ValueType>::EY;
            using WavefieldsEM<ValueType>::EZ;
            using WavefieldsEM<ValueType>::RX;
            using WavefieldsEM<ValueType>::RY;
            using WavefieldsEM<ValueType>::RZ;

            /* non-required wavefields */

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
