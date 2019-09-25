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
include(CMakePackageConfigHelpers)
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


### get git commit ###
execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
add_definitions("-DBUILD_COMMIT=${GIT_COMMIT}")


####################################################
#  Find LAMA library (mandatory)                   #
####################################################

set ( SCAI_PATH_HINTS "${SCAI_DIR}" "${SCAI_HOME}" "$ENV{SCAI_HOME}" )

if ( ${CMAKE_VERSION} VERSION_GREATER 3.11)
    ## use new policy to configure with SCAI_ROOT, suppress warning
    cmake_policy( SET CMP0074 NEW )
else ()
    ## simulate new policy by adding SCAI_ROOT to hints
    set ( SCAI_PATH_HINTS ${SCAI_PATH_HINTS} "${SCAI_ROOT}" "$ENV{SCAI_ROOT}" )
endif ()


find_package( SCAI REQUIRED CONFIG 
              NO_CMAKE_PATH 
              NO_SYSTEM_ENVIRONMENT_PATH 
              PATHS ${SCAI_PATH_HINTS} )

## take over compilfer flags as used by LAMA

include_directories( ${SCAI_INCLUDE_DIRS} )

set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SCAI_CXX_FLAGS} )
set( SOFI_used_libs ${SCAI_LIBRARIES} )


####################################################
#  Find Geographer library (optional)              #
####################################################

set ( Geographer_PATH_HINTS "${Geographer_DIR}" "${Geographer_HOME}" "$ENV{Geograhper_HOME}" )

if ( ${CMAKE_VERSION} VERSION_GREATER 3.11)
    ## use new policy to configure with SCAI_ROOT, suppress warning
    cmake_policy( SET CMP0074 NEW )
else ()
    ## simulate new policy by adding SCAI_ROOT to hints
    set ( Geograhper_PATH_HINTS ${Geographer_PATH_HINTS} "${Geographer_ROOT}" "$ENV{Geographer_ROOT}" )
endif ()

find_package( Geographer CONFIG
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH 
              PATHS ${Geographer_PATH_HINTS} )
message( "Geographer_LIBRARIES= \t" ${Geographer_LIBRARIES}$ )

set(CMAKE_VERBOSE_MAKEFILE ON)

#paths for supermuc; TODO: turn them to input arguments
FIND_PACKAGE( Zoltan2 PATHS "/home/hpc/pr87si/di36qin/Code/Trilinos/build/packages/zoltan2/")
set(METIS_BASE "/lrz/sys/libraries/metis/5.1.0/i8r8")
set(PARMETIS_BASE  "/lrz/sys/libraries/parmetis/4.0.3/ibmmpi/i8r8")
find_library( METIS_LIB metis "${METIS_BASE}/lib" )
find_library( PARMETIS_LIB parmetis "${PARMETIS_BASE}/lib" )

if ( ${Geographer_FOUND} )
    include_directories( ${Geographer_INCLUDE_DIRS} )
    add_definitions( "-DUSE_GEOGRAPHER" )
    set( SOFI_used_libs ${Geographer_LIBRARIES} ${SOFI_used_libs} )

	if( Zoltan2_FOUND AND METIS_BASE AND PARMETIS_BASE  )
		link_directories( ${Zoltan2_LIBRARY_DIRS} ${METIS_BASE}/lib ${PARMETIS_BASE}/lib )
		message( "Zoltan2_LIBRARY_DIRS= " ${Zoltan2_LIBRARY_DIRS} )
	    add_definitions( "-DUSE_GEOGRAPHER_WRAPPERS" )
		### include parmetis and metis ###
		#target_include_directories( prog PUBLIC "${METIS_BASE}/include" )
		#target_include_directories( prog PUBLIC "${PARMETIS_BASE}/include" )
	endif()
endif ()

####################################################
#  Add CXX Flags                                   #
####################################################

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror" )

###################################################
#  Collect source files of SOFI library            #
#  and copy header files to install                #
####################################################

set ( SOFI_LIB_DIRS 
        Acquisition 
        Configuration  
        ForwardSolver   
        Wavefields
        CheckParameter 
        IO
        Common
        Filter
        Modelparameter  
        Partitioning )

foreach ( dir ${SOFI_LIB_DIRS} ) 

    file( GLOB_RECURSE files ${dir}/*.cpp ${dir}/*.hpp )
    set ( SOFI_FILES ${SOFI_FILES} ${files} )

    file( GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR}/${dir} ${dir}/*.hpp )
    foreach ( file ${files} )
	    configure_file( ${dir}/${file} include/${dir}/${file} COPYONLY )
    endforeach ()

endforeach ()

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)


####################################################
#  SOFI library                                    #
####################################################

include_directories( Common Configuration Acquisition Modelparameter )

add_library( SOFI ${SOFI_FILES} )

target_link_libraries( SOFI ${SOFI_used_libs} )

install( TARGETS SOFI DESTINATION lib )

####################################################
#  SOFI executable                                 #
####################################################

add_executable( prog SOFI.cpp )

target_link_libraries( prog SOFI ${SOFI_used_libs} )

set_target_properties( prog PROPERTIES OUTPUT_NAME SOFI )

install( TARGETS prog DESTINATION bin )
install( DIRECTORY ${CMAKE_BINARY_DIR}/include DESTINATION . )

####################################################
#  Unit Test                                       #
####################################################

find_package( GTest )

if ( ${GTEST_FOUND} )

   file( GLOB SOFI_UNIT_TEST_SOURCES Tests/UnitTest/*.cpp )

   add_executable( utest ${SOFI_UNIT_TEST_SOURCES} )

   target_include_directories( utest PRIVATE ForwardSolver Configuration Acquisition Modelparameter ${GTEST_INCLUDE_DIR} )
   target_link_libraries( utest ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} SOFI ${SOFI_used_libs} )

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

target_link_libraries( itest SOFI ${SOFI_used_libs} )

set_target_properties( itest PROPERTIES OUTPUT_NAME Test_integration )

install( TARGETS itest DESTINATION bin )

####################################################
#  Create Model                                    #
####################################################

add_executable( model Tools/CreateModel/TwoLayer.cpp )

target_link_libraries( model SOFI ${SOFI_used_libs} )

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

### support for FindPackage ###

configure_file ( SOFIConfig.cmake.in ${CMAKE_BINARY_DIR}/SOFIConfig.cmake @ONLY)
install ( FILES ${CMAKE_BINARY_DIR}/SOFIConfig.cmake DESTINATION share/cmake/SOFI )

