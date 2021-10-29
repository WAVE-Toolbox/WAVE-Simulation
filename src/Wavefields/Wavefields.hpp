#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Vector.hpp>

#include "../Common/HostPrint.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/hmemo/HArray.hpp>

namespace KITGPI
{

    //! \brief Wavefields namespace
    namespace Wavefields
    {

        /*! \brief Abstract class to handle the wavefields for the forward modelling.
         *
         * Wavefields implements some methods, which are required by all derived classes.
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

            /* Common */
            //! \brief Declare Wavefield pointer
            typedef std::shared_ptr<Wavefields<ValueType>> WavefieldPtr;
            
            //! Reset wavefields
            virtual void resetWavefields() = 0;

            virtual int getNumDimension() const = 0;
            virtual std::string getEquationType() const = 0;

            //! Declare getter variable for context pointer
            virtual scai::hmemo::ContextPtr getContextPtr() = 0;

            //! \brief Initialization
            virtual void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) = 0;
            
            virtual bool isFinite(scai::dmemo::DistributionPtr dist) = 0;

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

            /* Seismic */
            virtual scai::lama::DenseVector<ValueType> &getRefVX() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVY() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVZ() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefP() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefPup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefPdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefPleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefPright() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVXup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVXdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVXleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVXright() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVYup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVYdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVYleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVYright() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVZup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVZdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVZleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefVZright() = 0;

            virtual scai::lama::DenseVector<ValueType> &getRefSxx() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefSyy() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefSzz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefSyz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefSxz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefSxy() = 0;

            virtual scai::lama::DenseVector<ValueType> &getRefRxx() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRyy() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRzz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRyz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRxz() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRxy() = 0;

            /* EM */
            virtual scai::lama::DenseVector<ValueType> &getRefHX() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefHY() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefHZ() = 0;

            virtual scai::lama::DenseVector<ValueType> &getRefEX() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEY() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEZ() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEXup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEYup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEZup() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEXdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEYdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEZdown() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEXleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEYleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEZleft() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEXright() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEYright() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefEZright() = 0;

            virtual scai::lama::DenseVector<ValueType> &getRefRX() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRY() = 0;
            virtual scai::lama::DenseVector<ValueType> &getRefRZ() = 0;
          protected:
            /* Common */
            void resetWavefield(scai::lama::DenseVector<ValueType> &vector);
            void initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

            typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
            SparseFormat transformMatrixYXZ;
            SparseFormat transformMatrixXZY;
            SparseFormat transformMatrixX;
            SparseFormat transformMatrixY;
            
            int numDimension;
            std::string equationType;

            /* Seismic */
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
        
            /* EM */
            scai::lama::DenseVector<ValueType> HX;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> HY;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> HZ;  //!< magnatic intensity Wavefield 
            scai::lama::DenseVector<ValueType> EX;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EY;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZ;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EXup;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EYup;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZup;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EXdown;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EYdown;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZdown;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EXleft;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EYleft;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZleft;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EXright;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EYright;  //!< electric intensity Wavefield
            scai::lama::DenseVector<ValueType> EZright;  //!< electric intensity Wavefield

            scai::lama::DenseVector<ValueType> RX; //!< Relaxation parameter memory varible
            scai::lama::DenseVector<ValueType> RY; //!< Relaxation parameter memory varible
            scai::lama::DenseVector<ValueType> RZ; //!< Relaxation parameter memory varible
        };
    }
}
