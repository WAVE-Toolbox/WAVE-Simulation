###########################################################################
#
#     Start the build of SOFI by the following command:
# 
#     SCAI_DIR=<path_to_lama> GTEST_ROOT=<path_to_googletest> cmake ../src 
#
#     make 
#     make install
#
#     default installation directory is the build directory
#    
###########################################################################

project( SOFI )

cmake_minimum_required ( VERSION 3.0 )

if ( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif ()

## Set default installation directory same as build directory

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    set ( CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

####################################################
#  Find LAMA library (mandatory)                   #
####################################################

## SCAI Lama library required, use SCAI_DIR or SCAI_ROOT (policy 0074) to find it

if ( ${CMAKE_VERSION} VERSION_GREATER 3.11)
    cmake_policy( SET CMP0074 NEW )
endif ()

find_package( SCAI REQUIRED CONFIG 
              NO_CMAKE_PATH 
              NO_SYSTEM_ENVIRONMENT_PATH )

## take over compilfer flags as used by LAMA

set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SCAI_CXX_FLAGS} )

####################################################
#  Collect source files of SOFI library            #
#  and copy header files to install                #
####################################################

set ( SOFI_LIB_DIRS 
        Acquisition 
        Configuration  
        ForwardSolver   
        PartitionedInOut 
        Wavefields
        CheckParameter 
        Common
        Filter
        Modelparameter  
        Partitioning )

foreach ( dir ${SOFI_LIB_DIRS} ) 

    file( GLOB_RECURSE files ${dir}/*.cpp ${dir}/*.hpp )
    set ( SOFI_FILES ${SOFI_FILES} ${files} )

    file( GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR}/${dir} ${dir}/*.hpp )
    foreach ( file ${files} )
        configure_file( ${dir}/${file} include/${file} COPYONLY )
    endforeach ()

endforeach ()

####################################################
#  SOFI library                                    #
####################################################

include_directories( ${SCAI_INCLUDE_DIRS} )

include_directories( Common Configuration Acquisition Modelparameter )

add_library( SOFI SHARED ${SOFI_FILES} )

target_link_libraries( SOFI ${SCAI_LIBRARIES} )

install( TARGETS SOFI DESTINATION lib )

####################################################
#  SOFI executable                                 #
####################################################

add_executable( prog SOFI.cpp )

target_link_libraries( prog SOFI )

set_target_properties( prog PROPERTIES OUTPUT_NAME SOFI )

install( TARGETS prog DESTINATION bin )

####################################################
#  Unit Test                                       #
####################################################

find_package( GTest )

if ( ${GTEST_FOUND} )

   file( GLOB SOFI_UNIT_TEST_SOURCES Tests/UnitTest/*.cpp )

   add_executable( utest ${SOFI_UNIT_TEST_SOURCES} )

   target_include_directories( utest PRIVATE ForwardSolver Configuration Acquisition Modelparameter ${GTEST_INCLUDE_DIR} )
   target_link_libraries( utest ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} SOFI )

   set_target_properties( utest PROPERTIES OUTPUT_NAME Test_unit )

   install( TARGETS utest DESTINATION bin )

else()

    message( WARNING "Skip configuration for utest, please set GTEST_ROOT to googletest installation directory" )

    add_custom_target (
        utest
        COMMAND echo "ATTENTION: cannot build utest as googletest not found, set GTEST_ROOT" 
    )

endif ()

####################################################
#  Integration Test                                #
####################################################

set ( SOFI_INTEGRATION_TEST_SOURCES
    Tests/IntegrationTest/Test_CompareSeismogram.cpp
)

add_executable( itest ${SOFI_INTEGRATION_TEST_SOURCES} )

target_link_libraries( itest SOFI )

set_target_properties( itest PROPERTIES OUTPUT_NAME Test_integration )

install( TARGETS itest DESTINATION bin )

####################################################
#  Create Model                                    #
####################################################

add_executable( model Tools/CreateModel/TwoLayer.cpp )

target_link_libraries( model SOFI )

set_target_properties( model PROPERTIES OUTPUT_NAME TwoLayer )

install( TARGETS model DESTINATION bin/tools )

####################################################
#  Doxygen documentation                           #
####################################################

find_package( Doxygen )

message ( STATUS "DOXYGEN_FOUND=${DOXYGEN_FOUND}, DOXYGEN_EXECUTABLE=${DOXYGEN_EXECUTABLE}" )

if ( DOXYGEN_FOUND )

    set ( DOXYGEN_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen" )

    configure_file ( "${CMAKE_SOURCE_DIR}/../doc/doxygen/Doxyfile" "${DOXYGEN_BINARY_DIR}/Doxyfile" @ONLY )

    file ( MAKE_DIRECTORY ${DOXYGEN_BINARY_DIR}/html )

    add_custom_command (
        OUTPUT ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_BINARY_DIR}/Doxyfile
        DEPENDS ${DOXYGEN_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${DOXYGEN_BINARY_DIR}
    )

    add_custom_target (
        doc
        DEPENDS ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMENT "Creating doxygen doc."
    )

    install ( DIRECTORY   ${DOXYGEN_BINARY_DIR}/html 
              DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc )

else ()

    add_custom_target (
        doc
        COMMAND echo "ATTENTION: doxygen not found, cannot build documentation" 
    )

endif ()

####################################################
#  Guide (latex/pdf) document                      #
####################################################

find_package ( LATEX COMPONENTS PDFLATEX BIBTEX )

message ( STATUS "LATEX_FOUND=${LATEX_FOUND}, PDFLATEX_COMPILER=${PDFLATEX_COMPILER}, BIBTEX_COMPILER=${BIBTEX_COMPILER}" )

if ( LATEX_FOUND )

    set ( GUIDE_BINARY_DIR "${CMAKE_BINARY_DIR}/doc/guide" )
    configure_file ( "${CMAKE_SOURCE_DIR}/../doc/guide/WAVE_guide_EN.tex" "${GUIDE_BINARY_DIR}/WAVE_guide_EN.tex" @ONLY )
    
    add_custom_command (
        OUTPUT ${GUIDE_BINARY_DIR}/WAVE_guide_EN.pdf
        COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
        COMMAND ${BIBTEX_COMPILER} WAVE_guide_EN
        COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
        COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
        DEPENDS ${GUIDE_BINARY_DIR}/WAVE_guide_EN.tex
        WORKING_DIRECTORY ${GUIDE_BINARY_DIR}
    )

    add_custom_target (
        pdf
        DEPENDS ${GUIDE_BINARY_DIR}/WAVE_guide_EN.pdf 
        COMMENT "Creating pdf guide."
    )

else ()

    add_custom_target (
        pdf
        COMMAND echo "ATTENTION: pdflatex/bibtex not found, cannot build user guide" 
    )

endif ()
