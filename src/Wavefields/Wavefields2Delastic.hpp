

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

#include "Wavefields.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"

namespace KITGPI
{

    namespace Wavefields
    {

        /*! \brief The class FD2Delastic holds the wavefields for 2D elastic simulation
         *
         */
        template <typename ValueType>
        class FD2Delastic : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD2Delastic(){};

            //! Default destructor
            ~FD2Delastic(){};

            explicit FD2Delastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefP() override;
            scai::lama::DenseVector<ValueType> &getRefRxx() override;
            scai::lama::DenseVector<ValueType> &getRefRyy() override;
            scai::lama::DenseVector<ValueType> &getRefRzz() override;
            scai::lama::DenseVector<ValueType> &getRefRyz() override;
            scai::lama::DenseVector<ValueType> &getRefRxz() override;
            scai::lama::DenseVector<ValueType> &getRefRxy() override;
            scai::lama::DenseVector<ValueType> &getRefVZ() override;
            scai::lama::DenseVector<ValueType> &getRefSzz() override;
            scai::lama::DenseVector<ValueType> &getRefSyz() override;
            scai::lama::DenseVector<ValueType> &getRefSxz() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Delastic<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Delastic<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Delastic<ValueType> operator*(KITGPI::Wavefields::FD2Delastic<ValueType> rhs);
            KITGPI::Wavefields::FD2Delastic<ValueType> operator*=(KITGPI::Wavefields::FD2Delastic<ValueType> rhs);

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType partitionedOut) override;
	    
	    void minusAssign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void timesAssign(ValueType rhs);
	    
          private:
	    void getCurl(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &curl, scai::lama::Vector<ValueType> const &SWaveModulus);
	    void getDiv(KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::lama::Vector<ValueType> &div, scai::lama::Vector<ValueType> const &SWaveModulus);
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::Sxx;
            using Wavefields<ValueType>::Syy;
            using Wavefields<ValueType>::Sxy;

            /* non-required wavefields */
            using Wavefields<ValueType>::P;
            using Wavefields<ValueType>::Szz;
            using Wavefields<ValueType>::Syz;
            using Wavefields<ValueType>::Sxz;
            using Wavefields<ValueType>::VZ;
            using Wavefields<ValueType>::Rxx;
            using Wavefields<ValueType>::Ryy;
            using Wavefields<ValueType>::Rzz;
            using Wavefields<ValueType>::Ryz;
            using Wavefields<ValueType>::Rxz;
            using Wavefields<ValueType>::Rxy;

            std::string type = "Elastic2D";
        };
    }
}
