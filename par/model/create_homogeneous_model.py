from writeMtxVector import write_mtx_vector
import numpy as np

## Input parameter
Model=1       # Model=0 (Acoustic), Model=1 (Elastic)
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
write_mtx_vector(filename+'.density.mtx',RHO)

if Model==0:
    Lambda=velocityP*velocityP*density #Lame acoustic
    LAMBDA=np.full(N,Lambda)
    write_mtx_vector(filename+'.lambda.mtx',LAMBDA)
elif Model==1:
    Lambda=(velocityP*velocityP - 2*velocityS*velocityS) *density #Lame acoustic
    Mu=(velocityS*velocityS)*density
    LAMBDA=np.full(N,Lambda)
    MU=np.full(N,Mu)
    write_mtx_vector(filename+'.vs.mtx',VS)
    write_mtx_vector(filename+'.lambda.mtx',LAMBDA)
    write_mtx_vector(filename+'.mu.mtx',MU)
else:
    print 'Please choose valid Model-Type.'
