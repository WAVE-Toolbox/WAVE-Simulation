## Preparatory steps

Since this finite-difference simulation code is based on the [LAMA framework](www.libama.org), you have to install this framework on your machine. A detailed documentation how to install LAMA is given on their website.

We successfully tested LAMA on different operating systems e.g. macOS (Sierra, El Capitan) and Linux (SUSE) as well as on a wide range of architectures e.g. HPC systems and GPUs.

Before the modelling code can be started, you have to set the paths to the installation of the [LAMA framework](www.libama.org):
- `export SCAI_ROOT=[PATH_TO_LAMA_BUILD]`

Moreover, you have to add the LAMA framework to the compiler search paths:
- `export DYLD_LIBRARY_PATH=${SCAI_ROOT}/lib:${DYLD_LIBRARY_PATH}`
- `export LD_LIBRARY_PATH=${SCAI_ROOT}/lib:${LD_LIBRARY_PATH}`

## Start the simulation
The simulation can be started by the example start file:

 ``source start_FDSimulation.sh``

Different kinds of parallelization are possible:

- Using OpenMP (shared-memory parallel):
 - Set `export OMP_NUM_THREADS=2`
- Using OpenMPI/IntelMPI (distributed-memory parallel):
 - Modify `start_FDSimulation.sh` to e.g. `mpirun -np 2 ./../bin/SOFI`

The standard configuration of `start_FDSimulation.sh` is using both kinds of parallelization.

## Run the tests
To test the proper functionality of the installation, you can run the build in unit and integration tests.
Since the [Google Test framework](https://github.com/google/googletest) is used for the unit tests, the environment variable `GTEST_DIR` has to be set to the location of the compiled Google Test library (`libgtest.*` and `libgtest_main.*`). The integration test has no special requirements apart the LAMA framework and MPI.

``source start_tests.sh``

If one of the test cases fails, the script will give an error message and will abort.

## Data processing for parallel I/O

Since the modelling code supports parallel I/O operations some pre- and post-processing steps maybe required to split or merge data. In order to use the binaries mentioned below, you have to build [LAMA](www.libama.org) with the option `BUILD_EXAMPLES=ON`.

#### Data pre-processing
- to partition a single vector-file `${SCAI_ROOT}/lama/examples/io/vectorRepartition.exe filename.mtx 1 filename.%r.mtx {NProcessors}`

#### Data post-processing
- Merge to a single vector-file `${SCAI_ROOT}/lama/examples/io/vectorRepartition.exe filename.%r.mtx {NProcessors} filename.mtx 1`
