#!/bin/bash
rm seismograms/seismogram.mtx
export OMP_NUM_THREADS=2
mpirun -np 2 ./../bin/SOFI3Dacoustic "configuration/configuration.txt"