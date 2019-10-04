#!/bin/bash

export SCAI_UNSUPPORTED=IGNORE

rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpiexec -n 28  ./../installation/bin/SOFI "configuration/mconfigurations/configuration_s305x305x305_vg1_gp1_fd0_geoKmeans_28.txt"
