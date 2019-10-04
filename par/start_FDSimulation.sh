#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpiexec -n 27 ./../bin/SOFI "configuration/configurationVarGrid.txt"
