#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Vector.hpp>

#include "../Common/HostPrint.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>
#include "../Wavefields/Wavefields.hpp"

namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        /*! \brief Abstract class to handle the wavefields for the forward modelling.
         *
         * WavefieldsEM implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class WavefieldsEM : public Wavefields<ValueType>
        {

          public:
            //! Default constructor
            WavefieldsEM(){};
            //! Default destructor
            ~WavefieldsEM(){};

            /* Common */
            //! Reset wavefields
            virtual void resetWavefields() = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            bool isFinite(scai::dmemo::DistributionPtr dist) override;
            
            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            virtual void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) = 0;

            //! Operator overloading
            virtual void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void timesAssign(ValueType rhs) = 0;
            
            virtual void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) = 0;
            
            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> &getRefVX() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVY() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVZ() override;
            virtual scai::lama::DenseVector<ValueType> &getRefP() override;
            virtual scai::lama::DenseVector<ValueType> &getRefPup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefPdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefPleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefPright() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVXup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVXdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVXleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVXright() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVYup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVYdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVYleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVYright() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVZup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVZdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVZleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefVZright() override;

            virtual scai::lama::DenseVector<ValueType> &getRefSxx() override;
            virtual scai::lama::DenseVector<ValueType> &getRefSyy() override;
            virtual scai::lama::DenseVector<ValueType> &getRefSzz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefSyz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefSxz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefSxy() override;

            virtual scai::lama::DenseVector<ValueType> &getRefRxx() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRyy() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRzz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRyz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRxz() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRxy() override;

            /* EM */
            virtual scai::lama::DenseVector<ValueType> &getRefHX() override;
            virtual scai::lama::DenseVector<ValueType> &getRefHY() override;
            virtual scai::lama::DenseVector<ValueType> &getRefHZ() override;

            virtual scai::lama::DenseVector<ValueType> &getRefEX() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEY() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEZ() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEXup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEYup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEZup() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEXdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEYdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEZdown() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEXleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEYleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEZleft() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEXright() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEYright() override;
            virtual scai::lama::DenseVector<ValueType> &getRefEZright() override;

            virtual scai::lama::DenseVector<ValueType> &getRefRX() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRY() override;
            virtual scai::lama::DenseVector<ValueType> &getRefRZ() override;
            
          protected:
              /* Common */            
            using Wavefields<ValueType>::transformMatrixYXZ;
            using Wavefields<ValueType>::transformMatrixXZY;
            using Wavefields<ValueType>::transformMatrixX;
            using Wavefields<ValueType>::transformMatrixY;

            using Wavefields<ValueType>::numDimension;
            using Wavefields<ValueType>::equationType;

            /* Seismic */
            using Wavefields<ValueType>::VX;  //!< Wavefield for velocity in x
            using Wavefields<ValueType>::VY;  //!< Wavefield for velocity in y
            using Wavefields<ValueType>::VZ;  //!< Wavefield for velocity in z
            using Wavefields<ValueType>::Sxx; //!< Wavefield
            using Wavefields<ValueType>::Syy; //!< Wavefield
            using Wavefields<ValueType>::Szz; //!< Wavefield
            using Wavefields<ValueType>::Syz; //!< Wavefield
            using Wavefields<ValueType>::Sxz; //!< Wavefield
            using Wavefields<ValueType>::Sxy; //!< Wavefield
            using Wavefields<ValueType>::P;   //!< Wavefield
            using Wavefields<ValueType>::Pup;   //!< Wavefield
            using Wavefields<ValueType>::Pdown;   //!< Wavefield
            using Wavefields<ValueType>::Pleft;   //!< Wavefield
            using Wavefields<ValueType>::Pright;   //!< Wavefield
            using Wavefields<ValueType>::VXup;
            using Wavefields<ValueType>::VXdown;
            using Wavefields<ValueType>::VXleft;
            using Wavefields<ValueType>::VXright;
            using Wavefields<ValueType>::VYup;
            using Wavefields<ValueType>::VYdown;
            using Wavefields<ValueType>::VYleft;
            using Wavefields<ValueType>::VYright;
            using Wavefields<ValueType>::VZup;
            using Wavefields<ValueType>::VZdown;
            using Wavefields<ValueType>::VZleft;
            using Wavefields<ValueType>::VZright;

            using Wavefields<ValueType>::Rxx; //!< Relaxation parameter
            using Wavefields<ValueType>::Ryy; //!< Relaxation parameter
            using Wavefields<ValueType>::Rzz; //!< Relaxation parameter
            using Wavefields<ValueType>::Ryz; //!< Relaxation parameter
            using Wavefields<ValueType>::Rxz; //!< Relaxation parameter
            using Wavefields<ValueType>::Rxy; //!< Relaxation parameter
            
            /* EM */
            using Wavefields<ValueType>::HX;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::HY;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::HZ;  //!< magnatic intensity Wavefield 
            using Wavefields<ValueType>::EX;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EY;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZ;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EXup;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EYup;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZup;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EXdown;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EYdown;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZdown;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EXleft;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EYleft;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZleft;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EXright;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EYright;  //!< electric intensity Wavefield
            using Wavefields<ValueType>::EZright;  //!< electric intensity Wavefield

            using Wavefields<ValueType>::RX; //!< Relaxation parameter memory varible
            using Wavefields<ValueType>::RY; //!< Relaxation parameter memory varible
            using Wavefields<ValueType>::RZ; //!< Relaxation parameter memory varible
        };
    }
}
