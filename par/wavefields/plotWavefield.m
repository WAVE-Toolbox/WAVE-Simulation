clearvars; close all;

addpath('../model/');
addpath('../configuration')
config=conf('../configuration/configuration.txt');

%% Define input parameter
component='P';
shotnr=1;

filename_base=config.getString('WavefieldFileName');
equationtype=config.getString('equationType');
dimension=config.getString('dimension');

%filename='wavefield.shot_0.Acoustic3D.P'; % File name of the model
filename=['../',filename_base,'.shot_',num2str(shotnr),'.',equationtype,dimension,'.',component]; % File name of the model
NX=config.getValue('NX');  % Number of grid points in X
NY=config.getValue('NY');  % Number of grid points in Y
NZ=config.getValue('NZ');  % Number of grid points in Z
DT=config.getValue('DT');
NTFirst=floor(config.getValue('tFirstSnapshot')/DT+0.5); % First Timestep
NTLast=floor(config.getValue('tLastSnapshot')/DT+0.5); %Last Timestep
NTint=floor(config.getValue('tIncSnapshot')/DT+0.5);  %Timestep Interval
DH=config.getValue('DH');   % Spatial grid sampling
tIncSnapshot=config.getValue('tIncSnapshot');

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


for ii=NTFirst:NTint:NTLast-NTint

model=readModelfromMtx([filename '.' num2str(ii) '.mtx'],NX,NY,NZ);
% model=readModelfromLMF([filename '.' num2str(ii) '.lmf'],NX,NY,NZ);

%% Plot
imagesc(model(:,:,LAYER))
colormap(seismic);
hold on
% plot(sourcex,sourcey,'m*')
% plot(receiverx,receivery,'gv')
caxis([-caxis_value caxis_value])
title( ['t = ' num2str(ii*DT) ' s'])
xlabel('X in gridpoints')
ylabel('Y in gridpoints')
%saveas(gcf,'vgWavefield.epsc')
pause(0.1)


end
