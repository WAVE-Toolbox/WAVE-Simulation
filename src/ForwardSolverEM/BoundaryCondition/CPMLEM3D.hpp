#pragma once

#include "../../ForwardSolver/BoundaryCondition/CPML.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Class for the usage of CPML boundaries for 3-D FD Simulations
            /*!
			 * Calculation of the CPML coefficients and the application of the cpml to the derivates for 3-D FD Simulations
			 *
			 */
            template <typename ValueType>
            class CPMLEM3D : public CPML<ValueType>
            {
              public:
                //! Default constructor
                CPMLEM3D(){};

                //! Default destructor
                ~CPMLEM3D(){};

                ValueType estimateMemory(scai::IndexType BoundaryWidth, scai::IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                void init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, scai::IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface);

                void resetCPML();

                void apply_ezx(scai::lama::DenseVector<ValueType> &ezx);
                void apply_eyx(scai::lama::DenseVector<ValueType> &eyx);
                void apply_ezy(scai::lama::DenseVector<ValueType> &ezy);
                void apply_exy(scai::lama::DenseVector<ValueType> &exy);
                void apply_eyz(scai::lama::DenseVector<ValueType> &eyz);
                void apply_exz(scai::lama::DenseVector<ValueType> &exz);
                void apply_hyx(scai::lama::DenseVector<ValueType> &hyx);
                void apply_hzx(scai::lama::DenseVector<ValueType> &hzx);
                void apply_hxy(scai::lama::DenseVector<ValueType> &hxy);
                void apply_hzy(scai::lama::DenseVector<ValueType> &hzy);
                void apply_hxz(scai::lama::DenseVector<ValueType> &hxz);
                void apply_hyz(scai::lama::DenseVector<ValueType> &hyz);

              private:
                using CPML<ValueType>::active;

                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_hyx; //!< CPML memory Variable
                VectorType psi_hzx; //!< CPML memory Variable
                VectorType psi_hxy; //!< CPML memory Variable
                VectorType psi_hzy; //!< CPML memory Variable
                VectorType psi_hxz; //!< CPML memory Variable
                VectorType psi_hyz; //!< CPML memory Variable

                VectorType psi_ezx; //!< CPML memory Variable
                VectorType psi_eyx; //!< CPML memory Variable
                VectorType psi_ezy; //!< CPML memory Variable
                VectorType psi_exy; //!< CPML memory Variable
                VectorType psi_eyz; //!< CPML memory Variable
                VectorType psi_exz; //!< CPML memory Variable

                VectorType a_x;      //!< CPML coefficient
                VectorType a_x_half; //!< CPML coefficient for staggered gridpoints
                VectorType a_y;      //!< CPML coefficient
                VectorType a_y_half; //!< CPML coefficient for staggered gridpoints
                VectorType a_z;      //!< CPML coefficient
                VectorType a_z_half; //!< CPML coefficient for staggered gridpoints

                VectorType b_x;      //!< CPML coefficient
                VectorType b_x_half; //!< CPML coefficient for staggered gridpoint
                VectorType b_y;      //!< CPML coefficient
                VectorType b_y_half; //!< CPML coefficient for staggered gridpoints
                VectorType b_z;      //!< CPML coefficient
                VectorType b_z_half; //!< CPML coefficient for staggered gridpoints
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
