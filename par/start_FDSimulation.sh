#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpirun -n 4 ./../build/bin/Simulation "configuration/configuration.txt"
