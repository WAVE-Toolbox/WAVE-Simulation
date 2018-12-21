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

            //! \brief 3-D acoustic free surface
            template <typename ValueType>
            class FreeSurfaceAcoustic : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceAcoustic(){};

                virtual ~FreeSurfaceAcoustic() = 0;

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ,Acquisition::Coordinates const &modelCoordinates, ValueType DT, ValueType DH) override;

              protected:
                using FreeSurface<ValueType>::active;

                scai::lama::DenseVector<ValueType> setSurfaceZero; //!< Vector, which sets the wavefields at the surface to zero
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
