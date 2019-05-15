#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE

rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpiexec -n 27  ./../build/bin/SOFI "configuration/configurationVarGrid.txt"
