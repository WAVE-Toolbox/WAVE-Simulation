#!/bin/bash
rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=2
mpirun-openmpi-mp -np 2 ./../bin/SOFI "configuration/configuration.txt"
