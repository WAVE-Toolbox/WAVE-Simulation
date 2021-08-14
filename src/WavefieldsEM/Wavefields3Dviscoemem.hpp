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
            using Wavefields<ValueType>::HX;
            using Wavefields<ValueType>::HY;
            using Wavefields<ValueType>::HZ;
            using Wavefields<ValueType>::EX;
            using Wavefields<ValueType>::EY;
            using Wavefields<ValueType>::EZ;
            using Wavefields<ValueType>::RX;
            using Wavefields<ValueType>::RY;
            using Wavefields<ValueType>::RZ;
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

            std::string type = EquationType+std::to_string(NumDimension)+"D";
        };
    }
}
