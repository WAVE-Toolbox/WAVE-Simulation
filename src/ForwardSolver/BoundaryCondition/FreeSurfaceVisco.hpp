#pragma once

#include "../../Common/HostPrint.hpp"
#include "../../Modelparameter/Modelparameter.hpp"
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

            //! \brief 3-D visco free surface
            template <typename ValueType>
            class FreeSurfaceVisco : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceVisco(){};

                //! Default destructor
                virtual ~FreeSurfaceVisco() = 0;

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ,Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, ValueType DH) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector<ValueType> &onePlusLtauP, scai::lama::Vector<ValueType> &onePlusLtauS, ValueType DT);

                void setMemoryVariableToZero(scai::lama::Vector<ValueType> &Ryy);

              protected:
                using FreeSurface<ValueType>::active;

                scai::lama::SparseVector<ValueType> setSurfaceZero;                  //!< Vector, which sets the wavefields at the surface to zero
                scai::lama::SparseVector<ValueType> selectHorizontalUpdate;          //!< Vector, which sets everything besides the free surface to zero
                scai::lama::SparseVector<ValueType> scaleStressVerticalUpdate;       //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                scai::lama::SparseVector<ValueType> scaleStressHorizontalUpdate;     //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                scai::lama::SparseVector<ValueType> scaleRelaxationVerticalUpdate;   //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                scai::lama::SparseVector<ValueType> scaleRelaxationHorizontalUpdate; //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivativesrelaxation mechanism
                scai::lama::SparseVector<ValueType> temp;                            //!< temporary Sparse Vector (lives as long as the forward solver)
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
