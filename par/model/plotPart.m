clearvars; close all;

%% Define input parameter
filename='../partitition.mtx'; % File name of the model
filenames= [ "../part_kmeans_1.mtx", "../part_kmeans_2.mtx", "../part_kmeans_3.mtx",  "../part_kmeans_4.mtx",  "../part_kmeans_5.mtx","../part_kmeans_6.mtx",  "../part_kmeans_7.mtx" ]

filename_w='../weights.mtx'; % File name of the model

NX=100;  % Number of grid points i0.0363n X
NY=100;  % Number of grid points in Y
NZ=1;  % Number of grid points in Zweig
DH=50;   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

weights=readModelfromMtx(filename_w,NX,NY,NZ);

%% Read model
% for f = 1:length(filenames)
%     part=readModelfromMtx(filenames(f),NX,NY,NZ);
%     
%     
%     %% Plot
%     figure
%     imagesc(X,Y,part(:,:,LAYER))
%     %set(gca,'yDir','normal'); 
%     colorbar
%     xlabel('X in meter')
%     ylabel('Y in meter')
% end

file_ml= [ "../mlRound_0.mtx", "../mlRound_1.mtx", "../mlRound_2.mtx", "../mlRound_3.mtx", "../mlRound_4.mtx", "../mlRound_5.mtx"]
for f = 1:1 %length(file_ml)
    part=readModelfromMtx(file_ml(f),NX,NY,NZ);
    
    figure
    imagesc(X,Y,part(:,:,LAYER))
    colorbar
    title("ml round"+f)
    xlabel('X in meter')
    ylabel('Y in meter')
end

figure
imagesc(X,Y,weights(:,:,LAYER))
colorbar
title("node weights")
xlabel('X in meter')
ylabel('Y in meter')

filenameGeom='../geomPart.mtx'
geomPart=readModelfromMtx(filenameGeom,NX,NY,NZ);
imagesc(X,Y,geomPart(:,:,LAYER))
colorbar
title('geomPart')
xlabel('X in meter')
ylabel('Y in meter')

figure
finalPart=readModelfromMtx(filename,NX,NY,NZ);
imagesc(X,Y,finalPart(:,:,LAYER))
colorbar
title("final partition")
xlabel('X in meter')
ylabel('Y in meter')

