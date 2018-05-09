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

            //! \brief Class for the usage of CPML boundaries for 3-D FD Simulations
            /*!
			 * Calculation of the CPML coefficients and the application of the cpml to the derivates for 3-D FD Simulations
			 *
			 */
            template <typename ValueType>
            class CPML3D : public CPML<ValueType>
            {
              public:
                //! Default constructor
                CPML3D(){};

                //! Default destructor
                ~CPML3D(){};

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DT, scai::IndexType DH, scai::IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface);

                void resetCPML();

                void apply_sxx_x(scai::lama::Vector<ValueType> &sxx_x);
                void apply_sxy_x(scai::lama::Vector<ValueType> &sxy_x);
                void apply_sxz_x(scai::lama::Vector<ValueType> &sxz_x);
                void apply_sxy_y(scai::lama::Vector<ValueType> &sxy_y);
                void apply_syy_y(scai::lama::Vector<ValueType> &syy_y);
                void apply_syz_y(scai::lama::Vector<ValueType> &syz_y);
                void apply_sxz_z(scai::lama::Vector<ValueType> &sxz_z);
                void apply_syz_z(scai::lama::Vector<ValueType> &syz_z);
                void apply_szz_z(scai::lama::Vector<ValueType> &szz_z);
                void apply_vxx(scai::lama::Vector<ValueType> &vxx);
                void apply_vyx(scai::lama::Vector<ValueType> &vyx);
                void apply_vzx(scai::lama::Vector<ValueType> &vzx);
                void apply_vxy(scai::lama::Vector<ValueType> &vxy);
                void apply_vyy(scai::lama::Vector<ValueType> &vyy);
                void apply_vzy(scai::lama::Vector<ValueType> &vzy);
                void apply_vxz(scai::lama::Vector<ValueType> &vxz);
                void apply_vyz(scai::lama::Vector<ValueType> &vyz);
                void apply_vzz(scai::lama::Vector<ValueType> &vzz);

              private:
		      
                using CPML<ValueType>::active;

                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_vxx; //!< CPML memory Variable
                VectorType psi_vyx; //!< CPML memory Variable
                VectorType psi_vzx; //!< CPML memory Variable
                VectorType psi_vxy; //!< CPML memory Variable
                VectorType psi_vyy; //!< CPML memory Variable
                VectorType psi_vzy; //!< CPML memory Variable
                VectorType psi_vxz; //!< CPML memory Variable
                VectorType psi_vyz; //!< CPML memory Variable
                VectorType psi_vzz; //!< CPML memory Variable

                VectorType psi_sxx_x; //!< CPML memory Variable
                VectorType psi_sxy_x; //!< CPML memory Variable
                VectorType psi_sxz_x; //!< CPML memory Variable
                VectorType psi_sxy_y; //!< CPML memory Variable
                VectorType psi_syy_y; //!< CPML memory Variable
                VectorType psi_syz_y; //!< CPML memory Variable
                VectorType psi_sxz_z; //!< CPML memory Variable
                VectorType psi_syz_z; //!< CPML memory Variable
                VectorType psi_szz_z; //!< CPML memory Variable

                VectorType k_x;      //!< CPML coefficient
                VectorType k_x_half; //!< CPML coefficient for staggered gridpoints
                VectorType k_y;      //!< CPML coefficient
                VectorType k_y_half; //!< CPML coefficient for staggered gridpoints
                VectorType k_z;      //!< CPML coefficient
                VectorType k_z_half; //!< CPML coefficient for staggered gridpoints

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
