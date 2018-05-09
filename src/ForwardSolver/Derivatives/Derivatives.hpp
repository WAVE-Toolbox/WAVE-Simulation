#pragma once
#include "../../Configuration/Configuration.hpp"
#include <scai/common/Stencil.hpp>
#include <scai/lama.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

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

            //            template <typename Type>
            //            class FreeSurface2Dsh;

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
                //                template <typename>
                //                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurface2Dsh;

                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceElastic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic;
                template <typename>
                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceVisco;
                //                template <typename>
                //                friend class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH;

                //! \brief Default constructor
                Derivatives() : spatialFDorder(0), useFreeSurface(false){};

                //! \brief Default destructor
                ~Derivatives(){};

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

                //! \brief Getter method for derivative matrix DyfFreeSurface
                virtual scai::lama::Matrix<ValueType> const &getDyfFreeSurface() const;
                //! \brief Getter method for derivative matrix DybFreeSurface
                virtual scai::lama::Matrix<ValueType> const &getDybFreeSurface() const;

                //! \brief Initialization
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) = 0;

                //! \brief Getter method for spatial FD-order
                scai::IndexType getSpatialFDorder() const;

              protected:
                //! \brief Initializsation of the derivative matrices
                /*!
                 *
                 \param dist Distribution of the wavefield
                 \param ctx Context
                 \param NX Total number of grid points in X
                 \param NY Total number of grid points in Y
                 \param NZ Total number of grid points in Z
                 \param DH Grid spacing (equidistant)
                 \param DT Temporal sampling interval
                 \param spatialFDorder FD-order of spatial stencils
                 \param comm Communicator
                 */
                virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH, ValueType DT, scai::IndexType spatialFDorder, scai::dmemo::CommunicatorPtr comm) = 0;

                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

                void setFDCoef(scai::IndexType spFDo);

                void calcDxf(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDyf(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDzf(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

                void calcDyb(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

                void calcDyfFreeSurface(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDybFreeSurface(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

                typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
                scai::lama::StencilMatrix<ValueType> Dxf;                    //!< Derivative matrix Dxf
                scai::lama::StencilMatrix<ValueType> Dyf;                    //!< Derivative matrix Dyf
                scai::lama::StencilMatrix<ValueType> Dzf;                    //!< Derivative matrix Dzf
                scai::lama::StencilMatrix<ValueType> Dxb;                    //!< Derivative matrix Dxb
                scai::lama::StencilMatrix<ValueType> Dyb;                    //!< Derivative matrix Dyb
                scai::lama::StencilMatrix<ValueType> Dzb;                    //!< Derivative matrix Dzb
                SparseFormat DyfFreeSurface;                                 //!< Derivative matrix DyfFreeSurface
                SparseFormat DybFreeSurface;                                 //!< Derivative matrix DybFreeSurface

                scai::IndexType spatialFDorder; //!< FD-Order of spatial derivative stencils

                bool useFreeSurface; //!< Switch to use free surface or not

                scai::hmemo::HArray<ValueType> FDCoef_f; //!< FD-coefficients forward
                scai::hmemo::HArray<ValueType> FDCoef_b; //!< FD-coefficients backward

                scai::common::Stencil1D<ValueType> stencilFD; // FD-stencil

              private:
                typedef void (Derivatives<ValueType>::*setRowElements_DPtr)(scai::IndexType, scai::IndexType &, scai::IndexType &, scai::hmemo::ReadAccess<ValueType> &, scai::hmemo::ReadAccess<ValueType> &, scai::hmemo::WriteAccess<scai::IndexType> &, scai::hmemo::WriteAccess<scai::IndexType> &, scai::hmemo::WriteAccess<ValueType> &, scai::IndexType, scai::IndexType, scai::IndexType); //!< Pointer to set elements functions

                typedef scai::IndexType (Derivatives<ValueType>::*calcNumberRowElements_DPtr)(scai::IndexType, scai::IndexType, scai::IndexType, scai::IndexType); //!< Pointer to counting elements functions

                void calcDerivativeMatrix(scai::lama::Matrix<ValueType> &D, calcNumberRowElements_DPtr calcNumberRowElements_D, setRowElements_DPtr setRowElements_D, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::DistributionPtr dist);

                scai::IndexType calcNumberRowElements_Dxf(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                scai::IndexType calcNumberRowElements_Dyf(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                scai::IndexType calcNumberRowElements_Dzf(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

                scai::IndexType calcNumberRowElements_Dyb(scai::IndexType rowNumber, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

                void setRowElements_Dxf(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                void setRowElements_Dyf(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                void setRowElements_Dzf(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                void setRowElements_Dyb(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

                void setRowElements_DyfFreeSurface(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
                void setRowElements_DybFreeSurface(scai::IndexType rowNumber, scai::IndexType &countJA, scai::IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<scai::IndexType> &csrJALocal, scai::hmemo::WriteAccess<scai::IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
