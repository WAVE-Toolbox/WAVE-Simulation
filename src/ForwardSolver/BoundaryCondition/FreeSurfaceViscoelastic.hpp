#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../Modelparameter/Modelparameter.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "../Derivatives/FDTD3D.hpp"
#include "FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D viscoelastic free surface
            template <typename ValueType>
            class FreeSurfaceViscoelastic : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceViscoelastic(){};

                //! Default destructor
                virtual ~FreeSurfaceViscoelastic(){};

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector<ValueType> &onePlusLtauP, scai::lama::Vector<ValueType> &onePlusLtauS, ValueType DT);

              protected:
                using FreeSurface<ValueType>::active;
                using FreeSurface<ValueType>::setZeroFreeSurface;
                scai::lama::SparseVector<ValueType> selectFreeSurface;               //!< Vector, which sets everything besides the free surface to zero
                scai::lama::SparseVector<ValueType> scaleStressVerticalUpdate;       //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                scai::lama::SparseVector<ValueType> scaleStressHorizontalUpdate;     //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                std::vector<scai::lama::SparseVector<ValueType>> scaleRelaxationVerticalUpdate;   //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                std::vector<scai::lama::SparseVector<ValueType>> scaleRelaxationHorizontalUpdate; //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivativesrelaxation mechanism
                scai::lama::SparseVector<ValueType> temp;                            //!< temporary Sparse Vector (lives as long as the forward solver)
                scai::IndexType numRelaxationMechanisms;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
