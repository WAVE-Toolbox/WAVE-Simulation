## Preparatory steps
- `export SCAI_ROOT=[PATH_TO_LAMA_BUILD]`
- `export DYLD_LIBRARY_PATH=${SCAI_ROOT}/lib`
- `export LD_LIBRARY_PATH=${SCAI_ROOT}/lib`

## Using OpenMP
- Set `export OMP_NUM_THREADS=1`

## Using OpenMPI/IntelMPI
- Modify `start_FDSimulation.sh` to e.g. `mpirun -np 4 ./../FDSimulation`

## Start the simulation
- `source start_FDSimulation.sh`
