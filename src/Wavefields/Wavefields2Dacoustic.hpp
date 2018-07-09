

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

        /*! \brief The class FD2Dacoustic holds the wavefields for 2D acoustic simulation
         *
         */
        template <typename ValueType>
        class FD2Dacoustic : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            FD2Dacoustic(){equationType="acoustic"; numDimension=2;};

            //! Default destructor
            ~FD2Dacoustic(){};

            explicit FD2Dacoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            void resetWavefields() override;

            int getNumDimension() const;
            std::string getEquationType() const;
            
            /* Getter routines for non-required wavefields: Will throw an error */
            scai::lama::DenseVector<ValueType> &getRefVZ() override;
            scai::lama::DenseVector<ValueType> &getRefSxx() override;
            scai::lama::DenseVector<ValueType> &getRefSyy() override;
            scai::lama::DenseVector<ValueType> &getRefSzz() override;
            scai::lama::DenseVector<ValueType> &getRefSyz() override;
            scai::lama::DenseVector<ValueType> &getRefSxz() override;
            scai::lama::DenseVector<ValueType> &getRefSxy() override;
            scai::lama::DenseVector<ValueType> &getRefRxx() override;
            scai::lama::DenseVector<ValueType> &getRefRyy() override;
            scai::lama::DenseVector<ValueType> &getRefRzz() override;
            scai::lama::DenseVector<ValueType> &getRefRyz() override;
            scai::lama::DenseVector<ValueType> &getRefRxz() override;
            scai::lama::DenseVector<ValueType> &getRefRxy() override;

            scai::hmemo::ContextPtr getContextPtr() override;

            /* Overloading Operators */
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(ValueType rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*=(ValueType rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs);
            KITGPI::Wavefields::FD2Dacoustic<ValueType> operator*=(KITGPI::Wavefields::FD2Dacoustic<ValueType> rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType partitionedOut) override;

            void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            void plusAssign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void assign(KITGPI::Wavefields::Wavefields<ValueType>  &rhs);
            void timesAssign(ValueType rhs);
    
          private:
              
            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType; 
              
            /* required wavefields */
            using Wavefields<ValueType>::VX;
            using Wavefields<ValueType>::VY;
            using Wavefields<ValueType>::P;

            /* non-required wavefields */
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

            std::string type = "Acoustic2D";
        };
    }
}
