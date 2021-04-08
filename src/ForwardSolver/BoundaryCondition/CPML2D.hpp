#pragma once

#include "CPML.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Class for the usage of CPML boundaries for 2-D FD Simulations
            /*!
			 * Calculation of the CPML coefficients and the application of the cpml to the derivates for 2-D FD Simulations
			 *
			 */
            template <typename ValueType>
            class CPML2D : public CPML<ValueType>
            {
              public:
                //! Default constructor
                CPML2D(){};

                //! Default destructor
                ~CPML2D(){};

                ValueType estimateMemory(scai::IndexType BoundaryWidth, scai::IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                void init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, scai::IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface);

                void resetCPML();

                void apply_sxx_x(scai::lama::DenseVector<ValueType> &sxx_x);
                void apply_sxy_x(scai::lama::DenseVector<ValueType> &sxy_x);
                void apply_sxz_x(scai::lama::DenseVector<ValueType> &sxz_x);
                void apply_sxy_y(scai::lama::DenseVector<ValueType> &sxy_y);
                void apply_syy_y(scai::lama::DenseVector<ValueType> &syy_y);
                void apply_syz_y(scai::lama::DenseVector<ValueType> &syz_y);
                void apply_vxx(scai::lama::DenseVector<ValueType> &vxx);
                void apply_vyx(scai::lama::DenseVector<ValueType> &vyx);
                void apply_vzx(scai::lama::DenseVector<ValueType> &vzx);
                void apply_vxy(scai::lama::DenseVector<ValueType> &vxy);
                void apply_vyy(scai::lama::DenseVector<ValueType> &vyy);
                void apply_vzy(scai::lama::DenseVector<ValueType> &vzy);

              private:
                // For the CPML Sparse Vectors and Dense Vectors can be declared. The code will run without any further changes.

                using CPML<ValueType>::active;
                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_vxx; //!< CPML memory Variable
                VectorType psi_vyx; //!< CPML memory Variable
                VectorType psi_vzx; //!< CPML memory Variable
                VectorType psi_vxy; //!< CPML memory Variable
                VectorType psi_vyy; //!< CPML memory Variable
                VectorType psi_vzy; //!< CPML memory Variable

                VectorType psi_sxx_x; //!< CPML memory Variable
                VectorType psi_sxy_x; //!< CPML memory Variable
                VectorType psi_sxz_x; //!< CPML memory Variable
                VectorType psi_sxy_y; //!< CPML memory Variable
                VectorType psi_syy_y; //!< CPML memory Variable
                VectorType psi_syz_y; //!< CPML memory Variable

                VectorType a_x;      //!< CPML coefficient
                VectorType a_y;      //!< CPML coefficient
                VectorType a_x_half; //!< CPML coefficient for staggered gridpoints
                VectorType a_y_half; //!< CPML coefficient for staggered gridpoints

                VectorType b_x;      //!< CPML coefficient
                VectorType b_y;      //!< CPML coefficient
                VectorType b_x_half; //!< CPML coefficient for staggered gridpoint
                VectorType b_y_half; //!< CPML coefficient for staggered gridpoints
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
