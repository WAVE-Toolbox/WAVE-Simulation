

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Vector.hpp>

#include "../Common/HostPrint.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>
#include "../Common/Hilbert.hpp"

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
            virtual scai::lama::DenseVector<ValueType> &getRefPup();
            virtual scai::lama::DenseVector<ValueType> &getRefPdown();
            virtual scai::lama::DenseVector<ValueType> &getRefPleft();
            virtual scai::lama::DenseVector<ValueType> &getRefPright();
            virtual scai::lama::DenseVector<ValueType> &getRefVXup();
            virtual scai::lama::DenseVector<ValueType> &getRefVXdown();
            virtual scai::lama::DenseVector<ValueType> &getRefVXleft();
            virtual scai::lama::DenseVector<ValueType> &getRefVXright();
            virtual scai::lama::DenseVector<ValueType> &getRefVYup();
            virtual scai::lama::DenseVector<ValueType> &getRefVYdown();
            virtual scai::lama::DenseVector<ValueType> &getRefVYleft();
            virtual scai::lama::DenseVector<ValueType> &getRefVYright();
            virtual scai::lama::DenseVector<ValueType> &getRefVZup();
            virtual scai::lama::DenseVector<ValueType> &getRefVZdown();
            virtual scai::lama::DenseVector<ValueType> &getRefVZleft();
            virtual scai::lama::DenseVector<ValueType> &getRefVZright();

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
            
            bool isFinite(scai::dmemo::DistributionPtr dist);

            ValueType getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numWavefields);

            //! \brief memory estimation
            virtual ValueType estimateMemory(scai::dmemo::DistributionPtr dist) = 0;

            ValueType getMemoryWavefield(scai::dmemo::DistributionPtr dist);

            virtual void write(scai::IndexType snapType, std::string baseName, scai::IndexType t, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, Modelparameter::Modelparameter<ValueType> const &model, scai::IndexType fileFormat) = 0;

            //! Operator overloading
            virtual void minusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void plusAssign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void assign(KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;
            virtual void timesAssign(ValueType rhs) = 0;
            
            virtual void applyTransform(scai::lama::CSRSparseMatrix<ValueType> lhs, KITGPI::Wavefields::Wavefields<ValueType> &rhs) = 0;        
            virtual void decompose(IndexType decomposeType, KITGPI::Wavefields::Wavefields<ValueType> &wavefieldsDerivative, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives) = 0;    
            
            void initTransformMatrixXYZ(IndexType decomposeType, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx);
            void calcTransformMatrixXYZ(KITGPI::Acquisition::Coordinates<ValueType> modelCoordinates);
            scai::lama::CSRSparseMatrix<ValueType> const &getTransformMatrixYXZ();
            scai::lama::CSRSparseMatrix<ValueType> const &getTransformMatrixXZY();

            KITGPI::Wavefields::Wavefields<ValueType> &operator=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator-=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator+=(KITGPI::Wavefields::Wavefields<ValueType> &rhs);
            KITGPI::Wavefields::Wavefields<ValueType> &operator*=(ValueType rhs);

          protected:
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
            SparseFormat transformMatrixYXZ;
            SparseFormat transformMatrixXZY;
            
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
            scai::lama::DenseVector<ValueType> Pup;   //!< Wavefield
            scai::lama::DenseVector<ValueType> Pdown;   //!< Wavefield
            scai::lama::DenseVector<ValueType> Pleft;   //!< Wavefield
            scai::lama::DenseVector<ValueType> Pright;   //!< Wavefield
            scai::lama::DenseVector<ValueType> VXup;
            scai::lama::DenseVector<ValueType> VXdown;
            scai::lama::DenseVector<ValueType> VXleft;
            scai::lama::DenseVector<ValueType> VXright;
            scai::lama::DenseVector<ValueType> VYup;
            scai::lama::DenseVector<ValueType> VYdown;
            scai::lama::DenseVector<ValueType> VYleft;
            scai::lama::DenseVector<ValueType> VYright;
            scai::lama::DenseVector<ValueType> VZup;
            scai::lama::DenseVector<ValueType> VZdown;
            scai::lama::DenseVector<ValueType> VZleft;
            scai::lama::DenseVector<ValueType> VZright;

            scai::lama::DenseVector<ValueType> Rxx; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryy; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rzz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Ryz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxz; //!< Relaxation parameter
            scai::lama::DenseVector<ValueType> Rxy; //!< Relaxation parameter
        };
    }
}
