## Preparatory steps

Since this finite-difference simulation code is based on the [LAMA framework](www.libama.org), you have to install this framework on your machine. A detailed documentation how to install LAMA is given on their website.

We successfully tested LAMA on different operating systems e.g. macOS (Sierra, El Capitan) and Linux (SUSE) as well as on a wide range of architectures e.g. HPC systems and GPUs.

Before the modelling code can be started, you have compile the code using cmake. Following example shows how to compile with 4 tasks
Further information can be found in `README.cmake` in the main directory

- ``mkdir build && cd build``
- ``SCAI_DIR="lama install directory" cmake ../src/ -DCMAKE_INSTALL_PREFIX=./``
- ``make install -j 4``


## Start the simulation
The simulation can be started by the example start file:

 ``source start_FDSimulation.sh``

Switch of Lama warnings and Tracing

`export SCAI_UNSUPPORTED=IGNORE`
` export SCAI_TRACE=OFF`

Different kinds of parallelization are possible:

- Using OpenMP (shared-memory parallel):
 - Set `export OMP_NUM_THREADS=2`
- Using OpenMPI/IntelMPI (distributed-memory parallel):
 - Modify `start_FDSimulation.sh` to e.g. `mpirun -np 2 ./../bin/SOFI`

The standard configuration of `start_FDSimulation.sh` is using both kinds of parallelization.

Run Code on GPUs
`export SCAI_ASYNCHRONOUS=2`
`export SCAI_DEVICE=0,1,2,3 (...)` (Comma seperated list of GPU devices on one node
`export SCAI_CONTEXT=CUDA`


## Run the tests
To test the proper functionality of the installation, you can run the build in unit and integration tests.
Since the [Google Test framework](https://github.com/google/googletest) is used for the unit tests, the environment variable `GTEST_DIR` has to be set to the location of the compiled Google Test library (`libgtest.*` and `libgtest_main.*`). The integration test has no special requirements apart the LAMA framework and MPI.

``source start_tests.sh``

If one of the test cases fails, the script will give an error message and will abort.

# Parallel IO
specify .lmf or .su as file type in the configuration file to use parallel in and output. 
