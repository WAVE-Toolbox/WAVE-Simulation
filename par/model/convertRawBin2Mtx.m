%% Converts a raw binary float matrix model to a mtx vector
% This script can be used to convert model that are saved as
% raw floating number in binary format to the mtx format.

% The major order used in WAVE-Simulation++ is X,Z,Y With Y=depth
% if your binary input model has a different major order. this script has
% to be modified 

% the result can be checked with the script plotModel.m

%% Clear memory
clearvars;
close all;

%% Input parameter
NX=100; % Number of grid points in X
NY=100; % Number of grid points in Y
NZ=100; % Number of grid points in Z (set to 1 for 2-D case)
filenameRawBin='model.bin'; % Filename input raw binary
filenameMtx='model.density.mtx'; % Filename output mtx

%% Convert
% Open raw binary file
[fid]=fopen(filenameRawBin,'r');

vector=fread(fid,NX*NY*NZ,'float');

%% uncomment and modify this part if your major order is different from X,Z,Y


% This example transforms a model with major order y,x,z to x,z,y

% matrix=reshape(vector,NY,NX,NZ);
% matrix=permute(matrix,[2 3 1]);
% vector=matrix(:);

%% Write matrix to mtx format
writeVector2mtx(filenameMtx,vector);

% the result can be checked with the script plotModel.m