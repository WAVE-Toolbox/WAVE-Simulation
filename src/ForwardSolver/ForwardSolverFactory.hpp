
#pragma once

#include "ForwardSolver.hpp"

#include "ForwardSolver2Dacoustic.hpp"
#include "ForwardSolver2Delastic.hpp"
#include "ForwardSolver2Dvisco.hpp"

#include "ForwardSolver3Dacoustic.hpp"
#include "ForwardSolver3Delastic.hpp"
#include "ForwardSolver3Dvisco.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief Factory class.
        template <typename ValueType>
        class Factory
        {
          public:
            //! \brief Declare ForwardSolver pointer
            typedef typename ForwardSolver<ValueType>::ForwardSolverPtr ForwardSolverPtr;

            Factory() = delete;
            Factory(Factory const &) = delete;
            void operator=(Factory const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ForwardSolverPtr Create(std::string dimension, std::string type);
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
