clc
close all
clear all

filename_m='./model.mtx'; % File name of the model
filename_x='./coordinatesX.mtx'; % File name of the x coordinates
filename_y='./coordinatesY.mtx'; % File name of the y coordinates
filename_z='./coordinatesZ.mtx'; % File name of the z coordinates

% z value of the 2D plane to display
layerZ=168;

%Ruduce data to show every "increment" point
increment=1;

%intrerpolate variable grid models
interpolate=1;

DH=8;

%% Read model
% 
x=readVectorfromMtx(filename_x);
y=readVectorfromMtx(filename_y);
z=readVectorfromMtx(filename_z);
model=readVectorfromMtx(filename_m);
  
if size(x) ~= size(model)
    error('Error size of coordinate vector differs from size of model')
end

x=x(1:increment:end);
y=y(1:increment:end);
z=z(1:increment:end);
model=model(1:increment:end);


indeces=find(z==layerZ);

model_layer=model(indeces);
x_layer=x(indeces);
y_layer=y(indeces);


if ~interpolate
    
    model_reg=zeros(max(y),max(x));
    for ii=1:length(y_layer)
       xi=x_layer(ii);
       yi=y_layer(ii);
    
       model_reg(yi+1,xi+1)=model_layer(ii);
    end
    
else

    [xq,yq] = meshgrid(0:1:max(x)-1, 0:1:max(y)-1);
    model_reg = griddata(x_layer,y_layer,model_layer,xq,yq);
    
end
 
x_reg=(0:max(x))*DH;
y_reg=(0:max(y))*DH;

figure
imagesc(x_reg,y_reg,model_reg)
title(filename_m)
xlabel('x in m')
ylabel('y in m')
