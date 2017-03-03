#pragma once

#include "../../Acquisition/Coordinates.hpp"
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

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, IndexType NX, IndexType NY, IndexType NZ, ValueType DT, ValueType DH) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model, scai::lama::Vector &onePlusLtauP, scai::lama::Vector &onePlusLtauS);

              protected:
                using FreeSurface<ValueType>::active;

                scai::lama::DenseVector<ValueType> setSurfaceZero;         //!< Vector, which sets the wavefields at the surface to zero
                scai::lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                scai::lama::DenseVector<ValueType> scaleStressHorizontalUpdate;     //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the stress update
                scai::lama::DenseVector<ValueType> scaleRelaxationHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter for the update of the relaxation mechanism
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
