#pragma once

#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Abstract class for the calculation and application of cpml boundaries
            template <typename ValueType>
            class CPML
            {
              public:
                //! \brief Default constructor
                CPML(){};

                //! \brief Default destructor
                ~CPML(){};

                //! init CPML coefficient vectors and CPML memory variables
                virtual void init(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, IndexType DH, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, bool useFreeSurface) = 0;

              protected:
                void resetVector(scai::lama::DenseVector<ValueType> &vector);

                void initVector(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);

                void SetCoeffCPML(scai::lama::DenseVector<ValueType> &a, scai::lama::DenseVector<ValueType> &b, scai::lama::DenseVector<ValueType> &kInv, scai::lama::DenseVector<ValueType> &a_half, scai::lama::DenseVector<ValueType> &b_half, scai::lama::DenseVector<ValueType> &kInv_half, IndexType coord,
                                  IndexType gdist, IndexType BoundaryWidth, ValueType NPower, ValueType KMaxCPML, ValueType CenterFrequencyCPML, ValueType VMaxCPML, IndexType i, ValueType DT, ValueType DH);

                void ResetCoeffFreeSurface(scai::lama::DenseVector<ValueType> &a, scai::lama::DenseVector<ValueType> &b, scai::lama::DenseVector<ValueType> &kInv,
                                           scai::lama::DenseVector<ValueType> &a_half, scai::lama::DenseVector<ValueType> &b_half, scai::lama::DenseVector<ValueType> &kInv_half,
                                           IndexType i);

                /*inline*/ void applyCPML(scai::lama::Vector &Vec, scai::lama::DenseVector<ValueType> &Psi, scai::lama::DenseVector<ValueType> &a, scai::lama::DenseVector<ValueType> &b, scai::lama::DenseVector<ValueType> &kInv);

                scai::lama::DenseVector<ValueType> psi_vxx; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vyx; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vzx; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vxy; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vyy; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vzy; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vxz; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vyz; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_vzz; //!< CPML memory Variable

                scai::lama::DenseVector<ValueType> psi_sxx_x; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_sxy_x; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_sxz_x; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_sxy_y; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_syy_y; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_syz_y; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_sxz_z; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_syz_z; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_szz_z; //!< CPML memory Variable

                scai::lama::DenseVector<ValueType> psi_p_x; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_p_y; //!< CPML memory Variable
                scai::lama::DenseVector<ValueType> psi_p_z; //!< CPML memory Variable

                scai::lama::DenseVector<ValueType> k_x; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> b_x; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> a_x; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> k_y; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> b_y; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> a_y; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> k_z; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> b_z; //!< CPML coefficient
                scai::lama::DenseVector<ValueType> a_z; //!< CPML coefficient

                scai::lama::DenseVector<ValueType> k_x_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> b_x_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> a_x_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> k_y_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> b_y_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> a_y_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> k_z_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> b_z_half; //!< CPML coefficient for staggered gridpoints
                scai::lama::DenseVector<ValueType> a_z_half; //!< CPML coefficient for staggered gridpoints

                scai::lama::DenseVector<ValueType> update_PmlTemp; //!< temporary vector for pml application
            };
        } /* end namespace BoundaryCondition  */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
