
#pragma once

#include "ForwardSolver.hpp"

#include "ForwardSolver2Dtmem.hpp"
#include "ForwardSolver2Demem.hpp"
#include "ForwardSolver2Dviscotmem.hpp"
#include "ForwardSolver2Dviscoemem.hpp"
#include "ForwardSolver3Demem.hpp"
#include "ForwardSolver3Dviscoemem.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief Factory class.
        template <typename ValueType>
        class FactoryEM
        {
          public:
            //! \brief Declare ForwardSolver pointer
            typedef typename ForwardSolverEM<ValueType>::ForwardSolverPtr ForwardSolverPtr;

            FactoryEM() = delete;
            FactoryEM(FactoryEM const &) = delete;
            void operator=(FactoryEM const &) = delete;

            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type ( emem)
             */
            static ForwardSolverPtr Create(std::string dimension, std::string type);
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
