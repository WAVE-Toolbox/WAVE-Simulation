clear all; close all; clc

addpath("../model");
addpath('../configuration')
config=conf('../configuration/configuration.txt');
%% Define input parameter
filename='wavefieldAcoustic2D.P'; % File name of the model
filename_x='../configuration/coordinatesX.mtx'; % File name of the model
filename_y='../configuration/coordinatesY.mtx'; % File name of the model
filename_z='../configuration/coordinatesZ.mtx'; % File name of the model

x=readVectorfromMtx(filename_x);
y=readVectorfromMtx(filename_y);
z=readVectorfromMtx(filename_z);

NX=config.getValue('NX');  % Number of grid points in X
NY=config.getValue('NY');  % Number of grid points in Y
NZ=config.getValue('NZ');  % Number of grid points in Z
NTFirst=floor(config.getValue('tFirstSnapshot')/DT+0.5); % First Timestep
NTLast=floor(config.getValue('tLastSnapshot')/DT+0.5); %Last Timestep
NTint=floor(config.getValue('tIncSnapshot')/DT+0.5);  %Timestep Interval
DH=config.getValue('DH');   % Spatial grid sampling
tIncSnapshot=config.getValue('tIncSnapshot');

% receiverx=[151 151];
% receivery=[103 202];
% sourcex=[151];
% sourcey=[73];




LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model

X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

figure('Position', [10 10 700 600])
caxis_value=1.0e-1;
load 'seismic.map'
colorbar


for ii=NTFirst:NTint:NTLast-NTint
%ii=500 

filenameii = [filename '.' num2str(ii) '.mtx'];
wavefield=readVectorfromMtx(filenameii);

 [xq,yq] = meshgrid(0:1:NX-1, 0:1:NY-1);
 vq = griddata(x,y,wavefield,xq,yq);

%% Plot
imagesc(vq)
colormap(seismic);
hold on
% plot(sourcex,sourcey,'m*')
% plot(receiverx,receivery,'gv')
caxis([-caxis_value caxis_value])
title( ['t = ' num2str(ii*DT) ' s'])
xlabel('X in gridpoints')
ylabel('Y in gridpoints')
%axis square
%saveas(gcf,'vgWavefield.epsc')
pause(0.8)


end