clearvars; close all;
addpath('../model/');

%% Define input parameter
%filename='wavefield.shot_0.Acoustic3D.P'; % File name of the model
filename='wavefield.shot_1.Acoustic3D.P'; % File name of the model
NX=133;  % Number of grid points in X
NY=112;  % Number of grid points in Y
NZ=105;  % Number of grid points in Z
NTFirst=0; % First Timestep
NTLast=950; %Last Timestep
NTint=50;  %Timestep Interval
DH=50;   % Spatial grid sampling
tIncSnapshot=0.002;

% receiverx=[151 151];
% receivery=[103 202];
% sourcex=[151];
% sourcey=[73];


LAYER=50; % Define layer of 3D model to display as 2D slice

%% Read model

X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

figure('Position', [200 200 700 600])
%figure
caxis_value=5.0e-2;
load 'seismic.map'
colorbar


for ii=NTFirst:NTint:NTLast


model=readModelfromMtx([filename '.' num2str(ii) '.mtx'],NX,NY,NZ);


%% Plot
imagesc(model(:,:,LAYER))
colormap(seismic);
hold on
% plot(sourcex,sourcey,'m*')
% plot(receiverx,receivery,'gv')
caxis([-caxis_value caxis_value])
title( ['t = ' num2str(ii*tIncSnapshot) ' s'])
xlabel('X in gridpoints')
ylabel('Y in gridpoints')
%saveas(gcf,'vgWavefield.epsc')
pause(0.1)


end