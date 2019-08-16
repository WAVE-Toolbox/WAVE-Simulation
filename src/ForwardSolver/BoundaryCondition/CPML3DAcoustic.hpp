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

                ValueType estimateMemory(scai::IndexType BoundaryWidth, scai::IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

                void init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, scai::IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface);

                void resetCPML();

                void apply_vxx(scai::lama::DenseVector<ValueType> &vxx);
                void apply_vyy(scai::lama::DenseVector<ValueType> &vyy);
                void apply_vzz(scai::lama::DenseVector<ValueType> &vzz);
                void apply_p_x(scai::lama::DenseVector<ValueType> &p_x);
                void apply_p_y(scai::lama::DenseVector<ValueType> &p_y);
                void apply_p_z(scai::lama::DenseVector<ValueType> &p_z);

              private:
                using CPML<ValueType>::active;
                typedef typename CPML<ValueType>::VectorType VectorType;

                VectorType psi_vxx; //!< CPML memory Variable
                VectorType psi_vyy; //!< CPML memory Variable
                VectorType psi_vzz; //!< CPML memory Variable

                VectorType psi_p_x; //!< CPML memory Variable
                VectorType psi_p_y; //!< CPML memory Variable
                VectorType psi_p_z; //!< CPML memory Variable

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
