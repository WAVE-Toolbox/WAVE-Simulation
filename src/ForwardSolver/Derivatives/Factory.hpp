

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
                static DerivativesPtr Create(std::string dimension)
                {

                    // transform to lower cases
                    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

                    // Assert correctness of input values
                    SCAI_ASSERT_ERROR(dimension.compare("3d") == 0 || dimension.compare("2d") == 0, "Unkown dimension");

                    if (dimension.compare("2d") == 0) {
                        return DerivativesPtr(new FDTD2D<ValueType>);
                    }
                    if (dimension.compare("3d") == 0) {
                        return DerivativesPtr(new FDTD3D<ValueType>);
                    }

                    COMMON_THROWEXCEPTION("Reached end of factory without match");
                    return nullptr;
                }
            };
        }
    }
}
