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

            explicit FD2Dviscoemem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefHX() override;
            scai::lama::DenseVector<ValueType> &getRefHY() override;
            scai::lama::DenseVector<ValueType> &getRefEZ() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRZ() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            ValueType estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs);
            KITGPI::Wavefields::FD2Dviscoemem<ValueType> operator*=(KITGPI::Wavefields::FD2Dviscoemem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            void timesAssign(scai::lama::DenseVector<ValueType> rhs);
            
            void applyTransform(scai::lama::Matrix<ValueType> const &lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposition, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            std::string EquationType;
            int NumDimension;
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;
            using Wavefields<ValueType>::numRelaxationMechanisms;

            /* required wavefields */
            using Wavefields<ValueType>::HZ;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::EX;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EY;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::RX;
            using Wavefields<ValueType>::RY;
            using Wavefields<ValueType>::EZup;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZdown;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZleft;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EZright;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EXup;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EXdown;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EXleft;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EXright;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EYup;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EYdown;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EYleft;  //!< electric intensity Wavefield 
            using Wavefields<ValueType>::EYright;  //!< electric intensity Wavefield

            /* non-required wavefields */
            using Wavefields<ValueType>::HX;
            using Wavefields<ValueType>::HY;
            using Wavefields<ValueType>::EZ;
            using Wavefields<ValueType>::RZ;

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
