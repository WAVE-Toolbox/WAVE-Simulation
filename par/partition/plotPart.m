clearvars; close all;


addpath("../model");

%% Define input parameter
filename='./partition.mtx'; % File name of the model
%filenames= [ "../part_kmeans_1.mtx", "../part_kmeans_2.mtx", "../part_kmeans_3.mtx",  "../part_kmeans_4.mtx",  "../part_kmeans_5.mtx","../part_kmeans_6.mtx",  "../part_kmeans_7.mtx" ]

filename_w='./weights.mtx'; % File name of the model
filename_x='../configuration/coordinatesX.mtx'; % File name of the model
filename_y='../configuration/coordinatesY.mtx'; % File name of the model
filename_z='../configuration/coordinatesZ.mtx'; % File name of the model

%Ruduce data to show every "increment" point
increment=1;


%% Read model

x=readVectorfromMtx(filename_x);
y=readVectorfromMtx(filename_y);
z=readVectorfromMtx(filename_z);

weights=readVectorfromMtx(filename_w);
part=readVectorfromMtx(filename);
  
x=x(1:increment:end);
y=y(1:increment:end);
z=z(1:increment:end);
part=part(1:increment:end);
weights=weights(1:increment:end);

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

figure 
histogram(part)

histo=zeros(1,max(part)+1);
for ii=0:max(part)
index=find(part==ii);
histo(ii+1)=sum(weights(index))/sum(weights);
end

figure
plot(0:max(part),histo)

