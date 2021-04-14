#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

rm -rf ../../WAVE-Inversion/par/seismograms/seismograms_checkerboard_dh10cm_Surface_StreamBig_Elastic2D_True.*.*.mtx
rm -rf wavefields/wavefield_checkerboard_dh10cm_Surface_StreamBig_True.shot_*.elastic2D.*.*.mtx
export OMP_NUM_THREADS=1
mpirun -np 6 ./../build/bin/Simulation "../../WAVE-Inversion/par/configuration/configuration_checkerboard_dh10cm_Surface_StreamBig_Elastic2D_True.txt"
