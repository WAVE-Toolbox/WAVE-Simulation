#pragma once

#include "../../Acquisition/Coordinates.hpp"
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

            //! \brief 3-D elastic free surface
            template <typename ValueType>
            class FreeSurfaceElastic : public FreeSurface<ValueType>
            {
              public:
                //! Default constructor
                FreeSurfaceElastic(){};

                virtual ~FreeSurfaceElastic() = 0;

                void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DT, ValueType DH) override;

                void setModelparameter(Modelparameter::Modelparameter<ValueType> const &model);

                void apply(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Syy, scai::lama::Vector<ValueType> &Szz);
                void apply(scai::lama::Vector<ValueType> &sumHorizonatlDerivative, scai::lama::Vector<ValueType> &Sxx, scai::lama::Vector<ValueType> &Syy);

              protected:
                using FreeSurface<ValueType>::active;

                scai::lama::DenseVector<ValueType> setSurfaceZero;         //!< Vector, which sets the wavefields at the surface to zero
                scai::lama::DenseVector<ValueType> selectHorizontalUpdate; //!< //!< Vector, which sets everything besides the free surface to zero

                scai::lama::DenseVector<ValueType> scaleHorizontalUpdate; //!< Vector, which sets the wavefields at the surface to zero which is scaled with the model parameter
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
