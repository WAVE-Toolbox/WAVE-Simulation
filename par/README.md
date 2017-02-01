## Preparatory steps

To start the modelling code, you have to set the paths to the LAMA framework:
- `export SCAI_ROOT=[PATH_TO_LAMA_BUILD]`
- `export DYLD_LIBRARY_PATH=${SCAI_ROOT}/lib:${DYLD_LIBRARY_PATH}`
- `export LD_LIBRARY_PATH=${SCAI_ROOT}/lib:${LD_LIBRARY_PATH}`

## Start the simulation
The simulation can be started by the example start file:

 ``source start_FDSimulation.sh``

Different kinds of parallelization are possible:

- Using OpenMP:
 - Set `export OMP_NUM_THREADS=1`
- Using OpenMPI/IntelMPI:
 - Modify `start_FDSimulation.sh` to e.g. `mpirun -np 4 ./../SOFI3Dacoustic`

## Start the tests
To test the proper functionality of the installation, you can run the tests:

- `source start_tests.sh`

## Data processing

Since the modelling code supports parallel I/O operations some pre- and post-processing steps maybe required to split or merge data.

#### Data pre-processing
- to partition a single vector-file `~/WAVE/scai_lama/build/lama/examples/io/vectorRepartition.exe filename.mtx 1 filename_%r.mtx {NProcessors}`

#### Data post-processing
- Merge to a single vector-file `~/WAVE/scai_lama/build/lama/examples/io/vectorRepartition.exe filename_%r.mtx {NProcessors} filename.mtx 1`
