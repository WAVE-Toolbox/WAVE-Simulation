#pragma once

#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurfaceVisco.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief 3-D visco free surface
            template <typename ValueType>
            class FreeSurface2Dvisco : public FreeSurfaceVisco<ValueType>
            {
              public:
                //! Default constructor
                FreeSurface2Dvisco(){};

                //! Default destructor
                ~FreeSurface2Dvisco(){};

                void apply(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &temp, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Syy, scai::lama::Vector<ValueType> &Rxx, scai::lama::Vector<ValueType> &Ryy);

              private:
                using FreeSurfaceVisco<ValueType>::setSurfaceZero;                  //!< Vector, which sets the wavefields at the surface to zero
                using FreeSurfaceVisco<ValueType>::scaleStressHorizontalUpdate;     //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the stress update
                using FreeSurfaceVisco<ValueType>::scaleRelaxationHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the update of the relaxation mechanism
                using FreeSurfaceVisco<ValueType>::active;
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
