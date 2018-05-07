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

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DT, scai::IndexType DH, scai::IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface);

                void resetCPML();

                void apply_sxx_x(scai::lama::Vector<ValueType> &sxx_x);
                void apply_sxy_x(scai::lama::Vector<ValueType> &sxy_x);
                void apply_sxy_y(scai::lama::Vector<ValueType> &sxy_y);
                void apply_syy_y(scai::lama::Vector<ValueType> &syy_y);
                void apply_vxx(scai::lama::Vector<ValueType> &vxx);
                void apply_vyx(scai::lama::Vector<ValueType> &vyx);
                void apply_vxy(scai::lama::Vector<ValueType> &vxy);
                void apply_vyy(scai::lama::Vector<ValueType> &vyy);

              private:
                // For the CPML Sparse Vectors and Dense Vectors can be declared. The code will run without any further changes.

                using CPML<ValueType>::active;
                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_vxx; //!< CPML memory Variable
                VectorType psi_vyx; //!< CPML memory Variable
                VectorType psi_vxy; //!< CPML memory Variable
                VectorType psi_vyy; //!< CPML memory Variable

                VectorType psi_sxx_x; //!< CPML memory Variable
                VectorType psi_sxy_x; //!< CPML memory Variable
                VectorType psi_sxy_y; //!< CPML memory Variable
                VectorType psi_syy_y; //!< CPML memory Variable

                VectorType k_x;      //!< CPML coefficient
                VectorType k_y;      //!< CPML coefficient
                VectorType k_x_half; //!< CPML coefficient for staggered gridpoints
                VectorType k_y_half; //!< CPML coefficient for staggered gridpoints

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
