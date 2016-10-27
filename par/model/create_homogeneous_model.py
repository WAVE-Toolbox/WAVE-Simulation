from writeMtxVector import write_mtx_vector
import numpy as np

## Input parameter
velocityP=3500 # P-wave velocity in m/s
velocityS=2000 # S-wave velocity in m/s
density=2000   # Density in Kg/m^3

N=100*100*100  # Number of grid points
filename='model' # Base filename

## Calculation

# Set homogeneous vectors
VP=np.full(N,velocityP)
VS=np.full(N,velocityS)
RHO=np.full(N,density)

write_mtx_vector(filename+'.vp.mtx',VP)
write_mtx_vector(filename+'.vs.mtx',VP)
write_mtx_vector(filename+'.density.mtx',VP)
