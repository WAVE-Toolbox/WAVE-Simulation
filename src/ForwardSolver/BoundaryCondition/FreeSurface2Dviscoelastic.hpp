#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurfaceViscoelastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D viscoelastic free surface
            template <typename ValueType>
            class FreeSurface2Dviscoelastic : public FreeSurfaceViscoelastic<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Dviscoelastic(){};

                //! Default destructor
                ~FreeSurface2Dviscoelastic(){};

                void exchangeHorizontalUpdate(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &vyy, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Rxx, ValueType DThalf);

              private:
                using FreeSurfaceViscoelastic<ValueType>::scaleStressHorizontalUpdate;     //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                using FreeSurfaceViscoelastic<ValueType>::scaleStressVerticalUpdate;       //!< Vector, which scales the horizontal stress updates to exchange vertical with horizontal derivatives
                using FreeSurfaceViscoelastic<ValueType>::scaleRelaxationHorizontalUpdate; //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                using FreeSurfaceViscoelastic<ValueType>::scaleRelaxationVerticalUpdate;   //!< Vector, which scales the horizontal relaxation updates to exchange vertical with horizontal derivatives
                using FreeSurfaceViscoelastic<ValueType>::selectFreeSurface;
                using FreeSurfaceViscoelastic<ValueType>::temp; //!< temporary Sparse Vector (lives as long as the forward solver)
                using FreeSurfaceViscoelastic<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
