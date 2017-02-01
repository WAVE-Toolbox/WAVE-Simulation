
#pragma once

#include "ForwardSolver.hpp"

#include "ForwardSolver2Dvisco.hpp"
#include "ForwardSolver2Delastic.hpp"
#include "ForwardSolver2Dacoustic.hpp"

#include "ForwardSolver3Dvisco.hpp"
#include "ForwardSolver3Delastic.hpp"
#include "ForwardSolver3Dacoustic.hpp"

namespace KITGPI {
    
    namespace ForwardSolver {
        
        //! \brief Factory class.
        template<typename ValueType>
        class Factory
        {
        public:
            
            //! \brief Declare ForwardSolver pointer
            typedef typename ForwardSolver<ValueType>::ForwardSolverPtr  ForwardSolverPtr;
            
            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            /*! \brief Create the right simmulation with factory methode.
             *
             \param dimension Dimension of the model (2D, 3D)
             \param type Simmulation type (acoustic, elsstic, viscoelastic)
             */
            static ForwardSolverPtr Create( std::string dimension, std::string type ) {
                
                // transform to lower cases
                std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
                std::transform(type.begin(), type.end(), type.begin(), ::tolower);
                
                // Assert correctness of input values
                SCAI_ASSERT_ERROR(dimension.compare("2d")==0 ||dimension.compare("3d")==0, "Unkown dimension" );
                SCAI_ASSERT_ERROR(type.compare("acoustic")==0 || type.compare("elastic")==0 || type.compare("visco")==0, "Unkown type" );
                
                // 2D
                if( dimension.compare("2d") == 0 && type.compare("acoustic") == 0 ){
                    return ForwardSolverPtr(new FD2Dacoustic<ValueType>);
                }
                if( dimension.compare("2d") == 0 && type.compare("elastic") == 0 ){
                    return ForwardSolverPtr(new FD2Delastic<ValueType>);
                }
                if( dimension.compare("2d") == 0 && type.compare("visco") == 0 ){
                    return ForwardSolverPtr(new FD2Dvisco<ValueType>);
                }
                
                // 3D
                if( dimension.compare("3d") == 0 && type.compare("acoustic") == 0 ){
                    return ForwardSolverPtr(new FD3Dacoustic<ValueType>);
                }
                if( dimension.compare("3d") == 0 && type.compare("elastic") == 0 ){
                    return ForwardSolverPtr(new FD3Delastic<ValueType>);
                }
                if( dimension.compare("3d") == 0 && type.compare("visco") == 0 ){
                    return ForwardSolverPtr(new FD3Dvisco<ValueType>);
                }
                
                COMMON_THROWEXCEPTION("Reached end of factory without match");
                return nullptr;
            }
            
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
