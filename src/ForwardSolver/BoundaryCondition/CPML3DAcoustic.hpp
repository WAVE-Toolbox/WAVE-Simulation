#pragma once

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
            class CPML3DAcoustic : public CPML<ValueType>
            {
              public:
                //! Default constructor
                CPML3DAcoustic(){};

                //! Default destructor
                ~CPML3DAcoustic(){};

                void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, scai::IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, scai::IndexType useFreeSurface);

                void resetCPML();

                void apply_vxx(scai::lama::Vector<ValueType> &vxx);
                void apply_vyy(scai::lama::Vector<ValueType> &vyy);
                void apply_vzz(scai::lama::Vector<ValueType> &vzz);
                void apply_p_x(scai::lama::Vector<ValueType> &p_x);
                void apply_p_y(scai::lama::Vector<ValueType> &p_y);
                void apply_p_z(scai::lama::Vector<ValueType> &p_z);

              private:
                using CPML<ValueType>::active;
                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_vxx; //!< CPML memory Variable
                VectorType psi_vyy; //!< CPML memory Variable
                VectorType psi_vzz; //!< CPML memory Variable

                VectorType psi_p_x; //!< CPML memory Variable
                VectorType psi_p_y; //!< CPML memory Variable
                VectorType psi_p_z; //!< CPML memory Variable

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

                //using CPML<ValueType>::applyCPML;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
