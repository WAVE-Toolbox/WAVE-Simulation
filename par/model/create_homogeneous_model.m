clearvars; close all;

%% Input parameter
velocityP=3500; % P-wave velocity in m/s
velocityS=2000; % S-wave velocity in m/s
density=2000; % Density in Kg/m^3

N=100*100*100; % Number of grid points
filename='model'; % Base filename

%% Calculation

% Set homogeneous vectors
VP=velocityP*ones(1,N);
VS=velocityS*ones(1,N);
RHO=density*ones(1,N);

% write to file
write_mtx_vector([filename '.vp.mtx'],VP);
write_mtx_vector([filename '.vs.mtx'],VS);
write_mtx_vector([filename '.density.mtx'],RHO);