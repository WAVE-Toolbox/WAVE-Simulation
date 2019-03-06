clearvars; close all;

%% Define input parameter
filename='../partitition.mtx'; % File name of the model
filename_w='../weights.mtx'; % File name of the model
filename_d='../damping.mtx'; % File name of the model
filename_x='../xcoords.mtx'; % File name of the model
filename_y='../ycoords.mtx'; % File name of the model
filename_z='../zcoords.mtx'; % File name of the model

NX=104;  % Number of grid points in X
NY=104;  % Number of grid points in Y
NZ=104;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model




part=readVectorfromMtx(filename);
weights=readVectorfromMtx(filename_w);
damping=readVectorfromMtx(filename_d);
x=readVectorfromMtx(filename_x);
y=readVectorfromMtx(filename_y);
z=readVectorfromMtx(filename_z);

% scatter3(x,y,z,50,damping,'filled','s')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% view(0,90)

% [xq,yq] = meshgrid(0:1:NX-1, 0:1:NY-1);
% vq = griddata(x,y,damping,xq,yq);
% 
% figure(1)
% imagesc(vq)
% axis ij

% [xq,yq,zq] = meshgrid(0:1:NX-1, 0:1:NY-1, 0:1:NZ-1);
% vq = griddata(x,y,z,part,xq,yq,zq);
% 
% 
% figure(1)
% imagesc(vq(:,:,60))
% %mesh(vq)
% axis ij

% X=0:DH:(NX*DH-DH);
% Y=0:DH:(NY*DH-DH);

% %% Plot
% figure
% imagesc(X,Y,part(:,:,LAYER))
% %set(gca,'yDir','normal'); 
% colorbar
% xlabel('X in meter')
% ylabel('Y in meter')
% 
% figure
% imagesc(X,Y,weights(:,:,LAYER))
% colorbar
% xlabel('X in meter')
% ylabel('Y in meter')