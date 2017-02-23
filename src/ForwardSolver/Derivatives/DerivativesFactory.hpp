

#pragma once

#include "Derivatives.hpp"
#include "FDTD2D.hpp"
#include "FDTD3D.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace Derivatives
        {

            //! \brief Factory class.
            template <typename ValueType>
            class Factory
            {

              public:
                //! \brief Declare Derivatives pointer
                typedef typename Derivatives<ValueType>::DerivativesPtr DerivativesPtr;

                Factory() = delete;
                Factory(Factory const &) = delete;
                void operator=(Factory const &) = delete;

                /*! \brief Create derivatives with factory methode.
                 *
                 \param dimension Dimension of the model (2D, 3D)
                 */
                static DerivativesPtr Create(std::string dimension);
            };
        }
    }
}
