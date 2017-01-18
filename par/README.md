## Preparatory steps
- `export SCAI_ROOT=[PATH_TO_LAMA_BUILD]`
- `export DYLD_LIBRARY_PATH=${SCAI_ROOT}/lib:${DYLD_LIBRARY_PATH}`
- `export LD_LIBRARY_PATH=${SCAI_ROOT}/lib:${LD_LIBRARY_PATH}`

## Using OpenMP
- Set `export OMP_NUM_THREADS=1`

## Using OpenMPI/IntelMPI
- Modify `start_FDSimulation.sh` to e.g. `mpirun -np 4 ./../SOFI3Dacoustic`

## Data pre-processing
- to partition a single vector-file `~/WAVE/scai_lama/build/lama/examples/io/vectorRepartition.exe filename.mtx 1 filename_%r.mtx {NProcessors}`

## Start the simulation
- `source start_FDSimulation.sh`

## Data post-processing
- Merge to a single vector-file `~/WAVE/scai_lama/build/lama/examples/io/vectorRepartition.exe filename_%r.mtx {NProcessors} filename.mtx 1`
