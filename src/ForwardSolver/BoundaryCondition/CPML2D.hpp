#pragma once

#include "../../Acquisition/Coordinates.hpp"
#include "../../Common/HostPrint.hpp"

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

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface);

                void reset();

                void apply_sxx_x(scai::lama::Vector &sxx_x);
                void apply_sxy_x(scai::lama::Vector &sxy_x);
                void apply_sxy_y(scai::lama::Vector &sxy_y);
                void apply_syy_y(scai::lama::Vector &syy_y);
                void apply_vxx(scai::lama::Vector &vxx);
                void apply_vyx(scai::lama::Vector &vyx);
                void apply_vxy(scai::lama::Vector &vxy);
                void apply_vyy(scai::lama::Vector &vyy);

              private:
                using CPML<ValueType>::psi_vxx;
                using CPML<ValueType>::psi_vyx;
                using CPML<ValueType>::psi_vxy;
                using CPML<ValueType>::psi_vyy;

                using CPML<ValueType>::psi_sxx_x;
                using CPML<ValueType>::psi_sxy_x;
                using CPML<ValueType>::psi_sxy_y;
                using CPML<ValueType>::psi_syy_y;

                using CPML<ValueType>::k_x;
                using CPML<ValueType>::b_x;
                using CPML<ValueType>::a_x;
                using CPML<ValueType>::k_y;
                using CPML<ValueType>::b_y;
                using CPML<ValueType>::a_y;

                using CPML<ValueType>::k_x_half;
                using CPML<ValueType>::b_x_half;
                using CPML<ValueType>::a_x_half;
                using CPML<ValueType>::k_y_half;
                using CPML<ValueType>::b_y_half;
                using CPML<ValueType>::a_y_half;

                /* non-required variables */
                using CPML<ValueType>::psi_vzx;
                using CPML<ValueType>::psi_vzy;
                using CPML<ValueType>::psi_vxz;
                using CPML<ValueType>::psi_vyz;
                using CPML<ValueType>::psi_vzz;
                using CPML<ValueType>::psi_sxz_x;
                using CPML<ValueType>::psi_syz_y;
                using CPML<ValueType>::psi_sxz_z;
                using CPML<ValueType>::psi_syz_z;
                using CPML<ValueType>::psi_szz_z;
                using CPML<ValueType>::k_z;
                using CPML<ValueType>::b_z;
                using CPML<ValueType>::a_z;
                using CPML<ValueType>::k_z_half;
                using CPML<ValueType>::b_z_half;
                using CPML<ValueType>::a_z_half;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
