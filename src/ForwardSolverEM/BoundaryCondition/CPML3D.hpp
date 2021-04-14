#pragma once

#include "../../Common/HostPrint.hpp"

#include "CPML.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Class for the usage of CPMLEM boundaries for 3-D FD Simulations
            /*!
			 * Calculation of the CPMLEM coefficients and the application of the cpml to the derivates for 3-D FD Simulations
			 *
			 */
            template <typename ValueType>
            class CPMLEM3D : public CPMLEM<ValueType>
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
                using CPMLEM<ValueType>::active;

                typedef typename CPMLEM<ValueType>::VectorType VectorType;

                VectorType psi_hyx; //!< CPMLEM memory Variable
                VectorType psi_hzx; //!< CPMLEM memory Variable
                VectorType psi_hxy; //!< CPMLEM memory Variable
                VectorType psi_hzy; //!< CPMLEM memory Variable
                VectorType psi_hxz; //!< CPMLEM memory Variable
                VectorType psi_hyz; //!< CPMLEM memory Variable

                VectorType psi_ezx; //!< CPMLEM memory Variable
                VectorType psi_eyx; //!< CPMLEM memory Variable
                VectorType psi_ezy; //!< CPMLEM memory Variable
                VectorType psi_exy; //!< CPMLEM memory Variable
                VectorType psi_eyz; //!< CPMLEM memory Variable
                VectorType psi_exz; //!< CPMLEM memory Variable

                VectorType a_x;      //!< CPMLEM coefficient
                VectorType a_x_half; //!< CPMLEM coefficient for staggered gridpoints
                VectorType a_y;      //!< CPMLEM coefficient
                VectorType a_y_half; //!< CPMLEM coefficient for staggered gridpoints
                VectorType a_z;      //!< CPMLEM coefficient
                VectorType a_z_half; //!< CPMLEM coefficient for staggered gridpoints

                VectorType b_x;      //!< CPMLEM coefficient
                VectorType b_x_half; //!< CPMLEM coefficient for staggered gridpoint
                VectorType b_y;      //!< CPMLEM coefficient
                VectorType b_y_half; //!< CPMLEM coefficient for staggered gridpoints
                VectorType b_z;      //!< CPMLEM coefficient
                VectorType b_z_half; //!< CPMLEM coefficient for staggered gridpoints
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
