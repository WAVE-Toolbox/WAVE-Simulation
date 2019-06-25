clearvars; close all;

%% Input parameter
velocityP=3500; % P-wave velocity in m/s
velocityS=2000; % S-wave velocity in m/s
density=2000; % Density in Kg/m^3

NX=100;
NY=100;
NZ=100; % Number of grid points
filename='model'; % Base filename

%% Calculation

% Set homogeneous vectors
VP=velocityP*ones(NY,NX,NZ);
VS=velocityS*ones(NY,NX,NZ);
RHO=density*ones(NY,NX,NZ);

% write to file
write3DModel2mtx([filename '.vp.mtx'],VP);
write3DModel2mtx([filename '.vs.mtx'],VS);
write3DModel2mtx([filename '.density.mtx'],RHO);