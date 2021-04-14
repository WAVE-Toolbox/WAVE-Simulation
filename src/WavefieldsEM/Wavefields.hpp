

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Vector.hpp>

#include "../Common/HostPrint.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../ModelparameterEM/Modelparameter.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

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
        class WavefieldsEM
        {

          public:
            //! Default constructor
            WavefieldsEM(){};
            //! Default destructor
            ~WavefieldsEM(){};

            //! \brief Declare Wavefield pointer
            typedef std::shared_ptr<WavefieldsEM<ValueType>> WavefieldPtr;

            //! Reset wavefields
            virtual void resetWavefields() = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            virtual scai::lama::DenseVector<ValueType> &getRefHX();
            virtual scai::lama::DenseVector<ValueType> &getRefHY();
            virtual scai::lama::DenseVector<ValueType> &getRefHZ();

            virtual scai::lama::DenseVector<ValueType> &getRefEX();
            virtual scai::lama::DenseVector<ValueType> &getRefEY();
            virtual scai::lama::DenseVector<ValueType> &getRefEZ();

            virtual scai::lama::DenseVector<ValueType> &getRefRX();
            virtual scai::lama::DenseVector<ValueType> &getRefRY();
            virtual scai::lama::DenseVector<ValueType> &getRefRZ();

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            bool isFinite(scai::dmemo::DistributionPtr dist);
            
            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numWavefields);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryWavefield(scai::dmemo::DistributionPtr dist);

            virtual void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::ModelparameterEM<ValueType> const &model, scai::IndexType fileFormat) = 0;

            //! Operator overloading
            virtual void minusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs) = 0;
            virtual void plusAssign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs) = 0;
            virtual void assign(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs) = 0;
            virtual void timesAssign(ValueType rhs) = 0;
            virtual void applyWavefieldTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs) = 0;

            KITGPI::Wavefields::WavefieldsEM<ValueType> &operator=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            KITGPI::Wavefields::WavefieldsEM<ValueType> &operator-=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            KITGPI::Wavefields::WavefieldsEM<ValueType> &operator+=(KITGPI::Wavefields::WavefieldsEM<ValueType> &rhs);
            KITGPI::Wavefields::WavefieldsEM<ValueType> &operator*=(ValueType rhs);

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            int numDimension;
            std::string equationType;

            scai::lama::DenseVector<ValueType> HX;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> HY;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> HZ;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> EX;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EY;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZ;  //!< electric intensity Wavefield

            scai::lama::DenseVector<ValueType> RX; //!< Relaxation parameter memory varible
            scai::lama::DenseVector<ValueType> RY; //!< Relaxation parameter memory varible
            scai::lama::DenseVector<ValueType> RZ; //!< Relaxation parameter memory varible
        };
    }
}
