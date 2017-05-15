#pragma once
#include "../../Configuration/Configuration.hpp"
#include <scai/lama.hpp>

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
                virtual scai::lama::Matrix const &getDxf() const;
                //! \brief Getter method for derivative matrix Dyf
                virtual scai::lama::Matrix const &getDyf() const;
                //! \brief Getter method for derivative matrix Dzf
                virtual scai::lama::Matrix const &getDzf() const;
                //! \brief Getter method for derivative matrix Dxb
                virtual scai::lama::Matrix const &getDxb() const;
                //! \brief Getter method for derivative matrix Dyb
                virtual scai::lama::Matrix const &getDyb() const;
                //! \brief Getter method for derivative matrix Dzb
                virtual scai::lama::Matrix const &getDzb() const;

                //! \brief Getter method for derivative matrix DyfVelocity
                virtual scai::lama::Matrix const &getDyfVelocity() const;
                //! \brief Getter method for derivative matrix DyfPressure
                virtual scai::lama::Matrix const &getDyfPressure() const;
                //! \brief Getter method for derivative matrix DybVelocity
                virtual scai::lama::Matrix const &getDybVelocity() const;
                //! \brief Getter method for derivative matrix DybPressure
                virtual scai::lama::Matrix const &getDybPressure() const;

                //! \brief Initialization
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm) = 0;

                //! \brief Getter method for spatial FD-order
                IndexType getSpatialFDorder() const;

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
                virtual void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, IndexType spatialFDorder, scai::dmemo::CommunicatorPtr comm) = 0;

                void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

                void setFDCoef(IndexType spFDo);

                void calcDxf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDyf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDzf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

                void calcDyb(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

                void calcDyfVelocity(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDyfPressure(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDybPressure(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);
                void calcDybVelocity(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

                typedef scai::lama::CSRSparseMatrix<ValueType> SparseFormat; //!< Define sparse format as CSRSparseMatrix
                SparseFormat Dxf;                                            //!< Derivative matrix Dxf
                SparseFormat Dyf;                                            //!< Derivative matrix Dyf
                SparseFormat Dzf;                                            //!< Derivative matrix Dzf
                SparseFormat Dxb;                                            //!< Derivative matrix Dxb
                SparseFormat Dyb;                                            //!< Derivative matrix Dyb
                SparseFormat Dzb;                                            //!< Derivative matrix Dzb
                SparseFormat DyfVelocity;                                    //!< Derivative matrix DyfVelocity
                SparseFormat DyfPressure;                                    //!< Derivative matrix DyfPressure
                SparseFormat DybVelocity;                                    //!< Derivative matrix DybVelocity
                SparseFormat DybPressure;                                    //!< Derivative matrix DybPressure

                IndexType spatialFDorder; //!< FD-Order of spatial derivative stencils

                bool useFreeSurface; //!< Switch to use free surface or not

                scai::hmemo::HArray<ValueType> FDCoef_f; //!< FD-coefficients forward
                scai::hmemo::HArray<ValueType> FDCoef_b; //!< FD-coefficients backward

              private:
                typedef void (Derivatives<ValueType>::*setRowElements_DPtr)(IndexType, IndexType &, IndexType &, scai::hmemo::ReadAccess<ValueType> &, scai::hmemo::ReadAccess<ValueType> &, scai::hmemo::WriteAccess<IndexType> &, scai::hmemo::WriteAccess<IndexType> &, scai::hmemo::WriteAccess<ValueType> &, IndexType, IndexType, IndexType); //!< Pointer to set elements functions

                typedef IndexType (Derivatives<ValueType>::*calcNumberRowElements_DPtr)(IndexType, IndexType, IndexType, IndexType); //!< Pointer to counting elements functions

                void calcDerivativeMatrix(scai::lama::Matrix &D, calcNumberRowElements_DPtr calcNumberRowElements_D, setRowElements_DPtr setRowElements_D, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist);

                IndexType calcNumberRowElements_Dxf(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
                IndexType calcNumberRowElements_Dyf(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);
                IndexType calcNumberRowElements_Dzf(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);

                IndexType calcNumberRowElements_Dyb(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ);

                void setRowElements_Dxf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
                void setRowElements_Dyf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
                void setRowElements_Dzf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
                void setRowElements_Dyb(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);

                void setRowElements_DyfVelocity(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
                void setRowElements_DybPressure(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
                void setRowElements_DybVelocity(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ);
            };
        } /* end namespace Derivatives */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
