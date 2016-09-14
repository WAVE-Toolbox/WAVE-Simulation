## Preparatory steps
- `export SCAI_ROOT=[PATH_TO_LAMA_BUILD]`
- `export DYLD_LIBRARY_PATH=${SCAI_ROOT}/lib`
- `export LD_LIBRARY_PATH=${SCAI_ROOT}/lib`

## Using OpenMP
- Set `export OMP_NUM_THREADS=1`

## Using OpenMPI/IntelMPI
- `mpirun -np 4 ./../FDSimulation`

## Run LAMA
- `source start_FDSimulation.sh`
