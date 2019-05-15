CMAKE and SOFI
==============

You can now use the CMake utility to build SOFI instead of using the Makefile as before.

Step 1: Configuration
---------------------

mkdir <build-directory>
cd <build-directory>
export SCAI_DIR=<lama-install-dir> 
export GTEST_ROOT=<googletest-install-dir>
cmake ../src

Note: use -DCMAKE_INSTALL_PREFIX=<install-dir> if you want to change the default installation directory
      (default is same as build directory)

# alternatively in one line
SCAI_DIR=<lama-installation-directory> GTEST_ROOT<googletest-installation-directory> cmake ../src 

# alternatively (with cmake version 3.12 or higher) use can use SCAI_ROOT instead of SCAI_DIR
cmake ../src -DSCAI_ROOT=<lama-install-dir> -DGTEST_ROOT=<googletest-install-dir>

Step 2: Build
-------------

make       ! generates library, SOFI executable, utest, itest
make doc   ! generates doxygen documentation
make pdf   ! generates pdf latex guide

Step 3: Install
---------------

make install    ! copies library to <install-dir>/lib, binaries to <install-dir>/bin

Note: Be default, installation directory is same as build directory.

Step 4: Run it

You find the executables in <install-dir>/bin  


Remarks:
========

- SCAI lama library is mandatory
- Googletest is optional, if not availabe or not found, the utest executable for unit tests is not built.
- the documentation will not be installed

Main differences to the old Makefile utility:
=============================================

- lib, utest, itest, model are always built as default
- generation of doxygen documentation and pdf guilde are separated
- make install must be called explicitly

Advantages of the new cmake utility:
====================================

- you have all the benefits that come with cmake

    make VERBOSE=1 ...    verbose information

- the file src/CMakeLists.txt is easier to understand than the old Makefile
- cmake can produce other build files than only for make
- separation between build and install makes a system-wide installation more convenient
- no platform dependencies
- cross-platform support
- finding 3rd party software is more convenient (e.g. doxygen, pdflatex, bibtex, googletest)

