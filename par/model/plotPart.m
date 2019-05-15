clearvars; close all;

%% Define input parameter
filename='../partitition.mtx'; % File name of the model
%filenames= [ "../part_kmeans_1.mtx", "../part_kmeans_2.mtx", "../part_kmeans_3.mtx",  "../part_kmeans_4.mtx",  "../part_kmeans_5.mtx","../part_kmeans_6.mtx",  "../part_kmeans_7.mtx" ]

filename_w='../weights.mtx'; % File name of the model
filename_d='../damping.mtx'; % File name of the model
filename_x='../configuration/coordinatesX.mtx'; % File name of the model
filename_y='../configuration/coordinatesY.mtx'; % File name of the model
filename_z='../configuration/coordinatesZ.mtx'; % File name of the model

filename_w='../weights.mtx'; % File name of the model
NX=100;  % Number of grid points in X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Z
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model

x=readVectorfromMtx(filename_x);
y=readVectorfromMtx(filename_y);
z=readVectorfromMtx(filename_z);

weights=readVectorfromMtx(filename_w);
part=readVectorfromMtx(filename);
  
figure  
scatter3(x,y,z,20,part,'filled','s')
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)


figure  
scatter3(x,y,z,20,weights,'filled','s')
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)

