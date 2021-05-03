#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF
export OMP_NUM_THREADS=1

rm -rf seismograms/seismogram.p.mtx

mpirun -n 4 ./../build/bin/Simulation "configuration/configuration.txt"
