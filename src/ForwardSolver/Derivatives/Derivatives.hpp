#pragma once
#include "../../Acquisition/Coordinates.hpp"
#include "../../Configuration/Configuration.hpp"
#include <map>
#include <scai/common/Stencil.hpp>
#include <scai/lama.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/lama/matrix/HybridMatrix.hpp>
#include <scai/tracing.hpp>

/* Forward declaration for friendship */
namespace KITGPI
{
    namespace ForwardSolver
    {
        namespace BoundaryCondition
        {
            template <typename Type>
            class FreeSurface2Delastic;

            template <typename Type>
            class FreeSurface2Dacoustic;

            template <typename Type>
            class FreeSurface2Dvisco;

            template <typename Type>
            class FreeSurface2Dsh;

            template <typename Type>
            class FreeSurfaceElastic;

            template <typename Type>
            class FreeSurfaceAcoustic;

            template <typename Type>
            class FreeSurfaceVisco;
        }
    }
}

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief Derivatives namespace
        namespace Derivatives
        {

            //! \brief Abstract class for the calculation of the Derivatives matrices
            template <typename ValueType>
            class Derivatives
            {
              public:
                //! \brief Declare Derivatives pointer
                typedef std::shared_ptr<Derivatives<ValueType>> DerivativesPtr;

                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Delastic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dacoustic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dvisco;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dsh;

                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco;

                //! \brief Default constructor
                Derivatives(): DyfFreeSurfaceHybrid(DyfFreeSurfaceStencil, DyfFreeSurfaceSparse),
                               DybFreeSurfaceHybrid(DybFreeSurfaceStencil, DybFreeSurfaceSparse){};

                //! \brief Default destructor
                ~Derivatives(){};

                void setup(Configuration::Configuration const &config);
                void setup(Configuration::Configuration const &config, std::vector<scai::IndexType> &FDorder);

                //! \brief Getter method for derivative matrix Dxf
                virtual scai::lama::Matrix<ValueType> const &getDxf() const;
                //! \brief Getter method for derivative matrix Dyf
                virtual scai::lama::Matrix<ValueType> const &getDyf() const;

                //! \brief Getter method for derivative matrix Dzf
                virtual scai::lama::Matrix<ValueType> const &getDzf() const;
                //! \brief Getter method for derivative matrix Dxb
                virtual scai::lama::Matrix<ValueType> const &getDxb() const;
                //! \brief Getter method for derivative matrix Dyb
                virtual scai::lama::Matrix<ValueType> const &getDyb() const;

                //! \brief Getter method for derivative matrix Dzb
                virtual scai::lama::Matrix<ValueType> const &getDzb() const;

                //! \brief Getter method for derivative matrix DyfStaggeredX
                virtual scai::lama::Matrix<ValueType> const &getDyfStaggeredX() const;
                //! \brief Getter method for derivative matrix DybStaggeredX
                virtual scai::lama::Matrix<ValueType> const &getDybStaggeredX() const;
                //! \brief Getter method for derivative matrix DyfStaggeredZ
                virtual scai::lama::Matrix<ValueType> const &getDyfStaggeredZ() const;
                //! \brief Getter method for derivative matrix DybStaggeredZ
                virtual scai::lama::Matrix<ValueType> const &getDybStaggeredZ() const;

                //! \brief Getter method for derivative matrix DyfFreeSurface
                virtual scai::lama::Matrix<ValueType> const &getDyfFreeSurface() const;
                scai::lama::Matrix<ValueType> &getDyfFreeSurface();
                
                //! \brief Getter method for derivative matrix DybFreeSurface
                virtual scai::lama::Matrix<ValueType> const &getDybFreeSurface() const;
                scai::lama::Matrix<ValueType> &getDybFreeSurface();
                
                virtual scai::lama::Matrix<ValueType> const &getDybStaggeredXFreeSurface() const;
                virtual scai::lama::Matrix<ValueType> const &getDybStaggeredZFreeSurface() const;

                //! \brief Getter method for interpolation matrix interFull
                virtual scai::lama::Matrix<ValueType> const *getInterFull() const;
                //! \brief Getter method for interpolation matrix interStaggeredX
                virtual scai::lama::Matrix<ValueType> const *getInterStaggeredX() const;
                //! \brief Getter method for interpolation matrix interStaggeredZ
                virtual scai::lama::Matrix<ValueType> const *getInterStaggeredZ() const;
                //! \brief Getter method for interpolation matrix interStaggeredXZ
                virtual scai::lama::Matrix<ValueType> const *getInterStaggeredXZ() const;

                virtual scai::lama::CSRSparseMatrix<ValueType> getCombinedMatrix() = 0;

                //! \brief Initialization
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) = 0;

                virtual void redistributeMatrices(scai::dmemo::DistributionPtr dist) = 0;

                virtual void estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

                //! \brief Getter method for spatial FD-order
                scai::IndexType getSpatialFDorder() const;

              protected:
                virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, ValueType DH, scai::dmemo::CommunicatorPtr comm) = 0;
                virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) = 0;
                virtual void initializeFreeSurfaceMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) = 0;

                void setFDCoef();

                void calcDxf(scai::dmemo::DistributionPtr dist);
                void calcDxf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDyf(scai::dmemo::DistributionPtr dist);
                void calcDyf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDzf(scai::dmemo::DistributionPtr dist);
                void calcDzf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDxb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDyb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDzb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);

                void calcDyfStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDybStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDyfStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDybStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);

                void calcDyfFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDybFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDybStaggeredXFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcDybStaggeredZFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);

                void calcInterpolationFull(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcInterpolationStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcInterpolationStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);
                void calcInterpolationStaggeredXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist);

                void setFDOrder(std::string const &FDorderFilename);

                void setFDOrder(std::vector<scai::IndexType> &FDorder);

                void setFDOrder(scai::IndexType FDorder);

                ValueType getMemoryStencilMatrix(scai::dmemo::DistributionPtr dist);
                ValueType getMemorySparseMatrix(scai::dmemo::DistributionPtr dist);
                ValueType getMemorySparseMatrix(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
                ValueType getMemoryInterpolationMatrix(scai::dmemo::DistributionPtr dist);

                typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
                scai::lama::StencilMatrix<ValueType> Dxf;                    //!< Derivative matrix Dxf
                scai::lama::StencilMatrix<ValueType> Dyf;                    //!< Derivative matrix Dyf
                scai::lama::StencilMatrix<ValueType> Dzf;                    //!< Derivative matrix Dzf
                scai::lama::StencilMatrix<ValueType> Dxb;                    //!< Derivative matrix Dxb
                scai::lama::StencilMatrix<ValueType> Dyb;                    //!< Derivative matrix Dyb
                scai::lama::StencilMatrix<ValueType> Dzb;                    //!< Derivative matrix Dzb

                SparseFormat DxfSparse; //!< Derivative matrix Dxf
                SparseFormat DyfSparse; //!< Derivative matrix Dyf
                SparseFormat DzfSparse; //!< Derivative matrix Dzf
                SparseFormat DxbSparse; //!< Derivative matrix Dxb
                SparseFormat DybSparse; //!< Derivative matrix Dyb
                SparseFormat DzbSparse; //!< Derivative matrix Dzb

                SparseFormat DyfStaggeredXSparse; //!< Derivative matrix Dyf for points staggered in x direction
                SparseFormat DybStaggeredXSparse; //!< Derivative matrix Dyf for points staggered in x direction
                SparseFormat DyfStaggeredZSparse; //!< Derivative matrix Dyf for points staggered in x direction
                SparseFormat DybStaggeredZSparse; //!< Derivative matrix Dyf for points staggered in x direction

//                 SparseFormat DyfFreeSurface; //!< Derivative matrix DyfFreeSurface
//                 SparseFormat DybFreeSurface; //!< Derivative matrix DybFreeSurface

                SparseFormat DyfFreeSurfaceSparse; //!< Derivative matrix DyfFreeSurface
                scai::lama::StencilMatrix<ValueType> DyfFreeSurfaceStencil;  //!< stencil part if hybrid format is used
                scai::lama::HybridMatrix<ValueType> DyfFreeSurfaceHybrid; // HybridMatrix( DyfFreeSurfaceSparse, DyfFreeSurfaceStencil );

                SparseFormat DybFreeSurfaceSparse; //!< Derivative matrix DybFreeSurface
                scai::lama::StencilMatrix<ValueType> DybFreeSurfaceStencil;  //!< stencil part if hybrid format is used
                scai::lama::HybridMatrix<ValueType> DybFreeSurfaceHybrid; // HybridMatrix( DybFreeSurfaceSparse, DybFreeSurfaceStencil );

                
                SparseFormat DybStaggeredXFreeSurface; //!< Derivative matrix DybFreeSurface
                SparseFormat DybStaggeredZFreeSurface; //!< Derivative matrix DybFreeSurface

                SparseFormat InterpolationFull;        //!< Interpolation matrix for P
                SparseFormat InterpolationStaggeredX;  //!< Interpolation matrix for staggerd points in x direction
                SparseFormat InterpolationStaggeredZ;  //!< Interpolation matrix for staggerd points in z direction
                SparseFormat InterpolationStaggeredXZ; //!< Interpolation matrix for staggerd points in x and in z direction

                ValueType DT;

                scai::IndexType useFreeSurface = 0; //!< Switch to use free surface or not
                bool useSparse = false;             //!< Switch to use Sparse Matrices
                bool useSparseFreeSurface = false; 
                bool useVarFDorder = false;         //!< Switch to use variable FDorder (layered)
                bool useVarGrid = false;            //!< Switch to use variable Grid
                bool isElastic = false;             //!< Switch to use variable Grid
              private:
                std::map<scai::IndexType, scai::common::Stencil1D<ValueType>> stencilFDmap; // FD-stencil
                                                                                            //     scai::IndexType spatialFDorder = 0;                                         //!< FD-Order of spatial derivative stencils
                std::vector<scai::IndexType> spatialFDorderVec;                             //!< std vector of variable FDordersof spatial derivative stencils  (layered)
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
