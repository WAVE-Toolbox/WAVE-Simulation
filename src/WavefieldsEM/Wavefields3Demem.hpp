

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

            explicit FD3Demem(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefRX() override;
            scai::lama::DenseVector<ValueType> &getRefRY() override;
            scai::lama::DenseVector<ValueType> &getRefRZ() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD3Demem<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*(KITGPI::Wavefields::FD3Demem<ValueType> rhs);
            KITGPI::Wavefields::FD3Demem<ValueType> operator*=(KITGPI::Wavefields::FD3Demem<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, scai::IndexType fileFormat) override;

            void minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            void timesAssign(ValueType rhs);
            
            void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs) override;
            void decompose(IndexType decomposeType, KITGPI::Wavefields::WavefieldsEM<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) override;

          private:
            void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &InverseDielectricPermittivity);
            void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &PWaveModulus);
            
            std::string EquationType;
            int NumDimension;
            using WavefieldsEM<ValueType>::numDimension;
            using WavefieldsEM<ValueType>::equationType;

            /* required wavefields */
            using WavefieldsEM<ValueType>::HX;  //!< magnatic intensity Wavefield 
            using WavefieldsEM<ValueType>::HY;  //!< magnatic intensity Wavefield 
            using WavefieldsEM<ValueType>::HZ;  //!< magnatic intensity Wavefield 
            using WavefieldsEM<ValueType>::EX;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EY;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EZ;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EZup;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EZdown;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EZleft;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EZright;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EXup;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EXdown;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EXleft;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EXright;  //!< electric intensity Wavefield
            using WavefieldsEM<ValueType>::EYup;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EYdown;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EYleft;  //!< electric intensity Wavefield 
            using WavefieldsEM<ValueType>::EYright;  //!< electric intensity Wavefield

            /* non-required wavefields */
            using WavefieldsEM<ValueType>::RX;
            using WavefieldsEM<ValueType>::RY;
            using WavefieldsEM<ValueType>::RZ;

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
