clearvars; %close all;

%% Define input parameter
filename='wavefieldAcoustic2D.VY'; % File name of the model
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
NT=950;
DH=50;   % Spatial grid sampling
tIncSnapshot=0.002;

receiver=[150 100 ; 150 200];
source=[150 70];


LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model

X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

figure('Position', [10 10 700 600])
caxis_value=5.0e-8;
%load 'seismic.map'
colorbar


for ii=0:50:NT;
%ii=500 

fileID = fopen([filename '.' num2str(ii) '.mtx'] ,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
A=fscanf(fileID,'%e',[1 size(1)*size(2)]);
model=reshape(A(:),[size(2), size(1)]);
model=permute(reshape(A,[NX,NY,NZ]),[2 1 3 4]);

max(max(model))

%% Plot
imagesc(model(:,:,LAYER))
%colormap(seismic);
caxis([-caxis_value caxis_value])
title( ['t = ' num2str(ii*tIncSnapshot) ' s'])
xlabel('X in gridpoints')
ylabel('Y in gridpoints')
axis square
saveas(gcf,'vgWavefield.epsc')
pause(0.1)

end