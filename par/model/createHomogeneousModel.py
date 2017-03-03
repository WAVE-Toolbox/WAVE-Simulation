from writeVector2mtx import writeVector2mtx
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

writeVector2mtx(filename+'.vp.mtx',VP)
writeVector2mtx(filename+'.density.mtx',RHO)

if Model==0:
    pWaveModulus=(velocityP*velocityP)*density #Lame acoustic
    PWAVEMODULUS=np.full(N,pWaveModulus)
    writeVector2mtx(filename+'.pWaveModulus.mtx',PWAVEMODULUS)
elif Model==1:
    pWaveModulus=(velocityP*velocityP)*density #Lame acoustic
    sWaveModulus=(velocityS*velocityS)*density
    PWAVEMODULUS=np.full(N,pWaveModulus)
    SWAVEMODULUS=np.full(N,sWaveModulus)
    writeVector2mtx(filename+'.vs.mtx',VS)
    writeVector2mtx(filename+'.pWaveModulus.mtx',PWAVEMODULUS)
    writeVector2mtx(filename+'.sWaveModulus.mtx',SWAVEMODULUS)
else:
    print 'Please choose valid Model-Type.'
