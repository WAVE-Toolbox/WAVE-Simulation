#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "WavefieldsEM.hpp"

namespace KITGPI
{

    namespace Wavefields
    {

        /*! \brief The class FD2Dtmem holds the wavefields for 2D tmem simulation
         *
         */
        template <typename ValueType>
        class FD2Dtmem : public WavefieldsEM<ValueType>
        {

          public:
            //! Default constructor
            FD2Dtmem():EquationType("tmem"),NumDimension(2)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD2Dtmem(){};

            explicit FD2Dtmem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefHZ() override;
            scai::lama::DenseVector<ValueType> &getRefEX() override;
            scai::lama::DenseVector<ValueType> &getRefEY() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRX() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRY() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRZ() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            ValueType estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dtmem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Dtmem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Dtmem<ValueType> operator*(KITGPI::Wavefields::FD2Dtmem<ValueType> rhs);
            KITGPI::Wavefields::FD2Dtmem<ValueType> operator*=(KITGPI::Wavefields::FD2Dtmem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void timesAssign(scai::lama::DenseVector<ValueType> rhs);
            
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity);

            std::string EquationType;
            int NumDimension;
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;
            using Wavefields<ValueType>::numRelaxationMechanisms;

            /* required wavefields */
            using Wavefields<ValueType>::HX;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::HY;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::EZ;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZup;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZdown;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZleft;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZright;  //!< electric intensity Wavefield

            /* non-required wavefields */
            using Wavefields<ValueType>::HZ;
            using Wavefields<ValueType>::EX;
            using Wavefields<ValueType>::EY;
            using Wavefields<ValueType>::RX;
            using Wavefields<ValueType>::RY;
            using Wavefields<ValueType>::RZ;

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
