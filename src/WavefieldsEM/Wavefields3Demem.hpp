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

        /*! \brief The class FD3Demem holds the wavefields for 3D emem simulation
         *
         */
        template <typename ValueType>
        class FD3Demem : public WavefieldsEM<ValueType>
        {

          public:
            //! Default constructor
            FD3Demem():EquationType("emem"),NumDimension(3)
            {
                equationType = EquationType;
                numDimension = NumDimension;
            };

            //! Default destructor
            ~FD3Demem(){};

            explicit FD3Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRX() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRY() override;
            std::vector<scai::lama::DenseVector<ValueType>> &getRefRZ() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            ValueType estimateMemory(dmemo::DistributionPtr dist, scai::IndexType numRelaxationMechanisms_in) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Demem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*(KITGPI::Wavefields::FD3Demem<ValueType> rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*=(KITGPI::Wavefields::FD3Demem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            
            void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) override;
            void decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &PWaveModulus);
            
            std::string EquationType;
            int NumDimension;
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;
            using Wavefields<ValueType>::numRelaxationMechanisms;

            /* required wavefields */
            using Wavefields<ValueType>::HX;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::HY;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::HZ;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::EX;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EY;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZ;  //!< electric intensity Wavefield
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
            using Wavefields<ValueType>::RX;
            using Wavefields<ValueType>::RY;
            using Wavefields<ValueType>::RZ;

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
