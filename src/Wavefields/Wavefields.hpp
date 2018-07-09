

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Vector.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>
#include "../Common/HostPrint.hpp"
#include "../PartitionedInOut/PartitionedInOut.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../Modelparameter/Modelparameter.hpp"


namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        /*! \brief Abstract class to handle the wavefields for the forward modelling.
         *
         * Wavefields implements some methods, which are requiered by all derived classes.
         * As this class is an abstract class, all methods are protected.
         */
        template <typename ValueType>
        class Wavefields
        {

          public:
            //! Default constructor
            Wavefields(){};
            //! Default destructor
            ~Wavefields(){};

            //! \brief Declare Wavefield pointer
            typedef std::shared_ptr<Wavefields<ValueType>> WavefieldPtr;

            //! Reset wavefields
            virtual void resetWavefields() = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;
            
            virtual scai::lama::DenseVector<ValueType> &getRefVX();
            virtual scai::lama::DenseVector<ValueType> &getRefVY();
            virtual scai::lama::DenseVector<ValueType> &getRefVZ();
            virtual scai::lama::DenseVector<ValueType> &getRefP();

            virtual scai::lama::DenseVector<ValueType> &getRefSxx();
            virtual scai::lama::DenseVector<ValueType> &getRefSyy();
            virtual scai::lama::DenseVector<ValueType> &getRefSzz();
            virtual scai::lama::DenseVector<ValueType> &getRefSyz();
            virtual scai::lama::DenseVector<ValueType> &getRefSxz();
            virtual scai::lama::DenseVector<ValueType> &getRefSxy();

            virtual scai::lama::DenseVector<ValueType> &getRefRxx();
            virtual scai::lama::DenseVector<ValueType> &getRefRyy();
            virtual scai::lama::DenseVector<ValueType> &getRefRzz();
            virtual scai::lama::DenseVector<ValueType> &getRefRyz();
            virtual scai::lama::DenseVector<ValueType> &getRefRxz();
            virtual scai::lama::DenseVector<ValueType> &getRefRxy();

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;

            virtual void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType partitionedOut) = 0;

            //! Operator overloading
            virtual void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void timesAssign(ValueType rhs) = 0;

            KITGPI::Wavefields::Wavefields<ValueType> &operator=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator-=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator+=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator*=(ValueType rhs);
	    
          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            void writeWavefield(scai::lama::Vector<ValueType> &vector, std::string component, std::string fileBaseName, scai::IndexType t, scai::IndexType partitionedOut);

            int numDimension;
            std::string equationType; 
            
            scai::lama::DenseVector<ValueType> VX;  //!< Wavefield for velocity in x
            scai::lama::DenseVector<ValueType> VY;  //!< Wavefield for velocity in y
            scai::lama::DenseVector<ValueType> VZ;  //!< Wavefield for velocity in z
            scai::lama::DenseVector<ValueType> Sxx; //!< Wavefield
            scai::lama::DenseVector<ValueType> Syy; //!< Wavefield
            scai::lama::DenseVector<ValueType> Szz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Syz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Sxz; //!< Wavefield
            scai::lama::DenseVector<ValueType> Sxy; //!< Wavefield
            scai::lama::DenseVector<ValueType> P;   //!< Wavefield

            scai::lama::DenseVector<ValueType> Rxx; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryy; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rzz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxy; //!< Relaxation parameter
        };
    }
}
