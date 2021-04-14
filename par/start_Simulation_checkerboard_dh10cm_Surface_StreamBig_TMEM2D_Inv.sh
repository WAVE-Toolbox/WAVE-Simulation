#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

rm -rf ../../WAVE-Inversion/par/seismograms/seismograms_checkerboard_dh10cm_Surface_StreamBig_TMEM2D_Start.*.*.mtx
rm -rf wavefields/wavefield_checkerboard_dh10cm_Surface_StreamBig_Start.shot_*.tmem2D.*.*.mtx
export OMP_NUM_THREADS=1
mpirun -np 6 ./../build/bin/Simulation "../../WAVE-Inversion/par/configuration/configuration_checkerboard_dh10cm_Surface_StreamBig_TMEM2D_Inv_HPC.txt"
