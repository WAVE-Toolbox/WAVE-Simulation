#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE

rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpirun -np 4  ./../build/bin/Simulation "configuration/configurationVarGrid.txt"
