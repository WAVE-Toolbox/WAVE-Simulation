#!/bin/bash
rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=2
mpirun -np 2 ./../build/bin/SOFI "configuration/configuration.txt"
