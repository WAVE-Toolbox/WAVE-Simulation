clearvars; close all;

%% Define input parameter
filename='model.vp.mtx'; % File name of the model
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=100;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model
model=readModelfromMtx(filename,NX,NY,NZ);
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
imagesc(X,Y,model(:,:,LAYER))
colorbar
xlabel('X in meter')
ylabel('Y in meter')