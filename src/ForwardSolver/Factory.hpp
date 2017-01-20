
#pragma once

#include "ForwardSolver.hpp"

#include "ForwardSolver2Dvisco.hpp"
#include "ForwardSolver2Delastic.hpp"
#include "ForwardSolver2Dacoustic.hpp"

#include "ForwardSolver3Dvisco.hpp"
#include "ForwardSolver3Delastic.hpp"
#include "ForwardSolver3Dacoustic.hpp"

namespace KITGPI {
    
    //! \brief ForwardSolver namespace
    namespace ForwardSolver {
        
        //! \brief Abstract class for forward solver
        template<typename ValueType>
        class Factory
        {
        public:
            
            typedef typename ForwardSolver<ValueType>::ForwardSolverPtr  ForwardSolverPtr;
            
            Factory() = delete;
            Factory(Factory const&) = delete;
            void operator=(Factory const&) = delete;
            
            static ForwardSolverPtr Create( std::string const& dimension, std::string const& type ) {
                
                // Assert correctness of input values
                SCAI_ASSERT_ERROR(dimension.compare("2D")==0 ||dimension.compare("3D")==0, "Unkown dimension" );
                SCAI_ASSERT_ERROR(type.compare("acoustic")==0 || type.compare("elastic")==0 || type.compare("visco")==0, "Unkown type" );
                
                // 2D
                if( dimension.compare("2D") == 0 && type.compare("acoustic") == 0 ){
                    return ForwardSolverPtr(new FD2Dacoustic<ValueType>);
                }
                if( dimension.compare("2D") == 0 && type.compare("elastic") == 0 ){
                    return ForwardSolverPtr(new FD2Delastic<ValueType>);
                }
                if( dimension.compare("2D") == 0 && type.compare("visco") == 0 ){
                    return ForwardSolverPtr(new FD2Dvisco<ValueType>);
                }
                
                // 3D
                if( dimension.compare("3D") == 0 && type.compare("acoustic") == 0 ){
                    return ForwardSolverPtr(new FD3Dacoustic<ValueType>);
                }
                if( dimension.compare("3D") == 0 && type.compare("elastic") == 0 ){
                    return ForwardSolverPtr(new FD3Delastic<ValueType>);
                }
                if( dimension.compare("3D") == 0 && type.compare("visco") == 0 ){
                    return ForwardSolverPtr(new FD3Dvisco<ValueType>);
                }
                
                COMMON_THROWEXCEPTION("Reached end of factory without match");
                return nullptr;
            }
            
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
