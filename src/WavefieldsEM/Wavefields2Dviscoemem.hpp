

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

        /*! \brief The class FD2Dviscoemem holds the wavefields for 2D viscoemem emem simulation
         *
         */
        template <typename ValueType>
        class FD2Dviscoemem : public WavefieldsEM<ValueType>
        {

          public:
            //! Default constructor
            FD2Dviscoemem():EquationType("viscoemem"),NumDimension(2)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD2Dviscoemem(){};

            explicit FD2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefHX() override;
            scai::lama::DenseVector<ValueType> &getRefHY() override;
            scai::lama::DenseVector<ValueType> &getRefEZ() override;
            scai::lama::DenseVector<ValueType> &getRefRZ() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*=(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void applyWavefieldTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);

          private:
            std::string EquationType;
            int NumDimension;
            using WavefieldsEM<ValueType>::numDimension;
            using WavefieldsEM<ValueType>::equationType;

            /* required wavefields */
            using WavefieldsEM<ValueType>::HZ;  //!< magnatic intensity Wavefield 
            using WavefieldsEM<ValueType>::EX;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EY;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::RX;
            using WavefieldsEM<ValueType>::RY;

            /* non-required wavefields */
            using WavefieldsEM<ValueType>::HX;
            using WavefieldsEM<ValueType>::HY;
            using WavefieldsEM<ValueType>::EZ;
            using WavefieldsEM<ValueType>::RZ;

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
