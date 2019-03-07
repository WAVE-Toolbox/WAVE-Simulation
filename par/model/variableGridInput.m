close all;
clear all
clc



%% Input
% Gridpoints of layer with dh=1 (every gridpoint)
NXmax=305;
NZmax=1;
NYmax=305;

interface(1)=100;
% interface(2)=151;
% interface(3)=202;


% interface(1)=10;
% interface(2)=37;
% interface(3)=55;

% 
% 
 dh(1)=1;
 dh(2)=3;
%  dh(3)=9;
%  dh(4)=1;


%% constructor
interface(length(interface)+1)=NYmax-1;


DHmax=max(dh);
if (NXmax ~= floor(NXmax/DHmax)*DHmax+1 + floor(DHmax/2))
    message ="NXmax must be " + num2str(floor(NXmax/DHmax)*DHmax+1 + floor(DHmax/2));
    disp(message)
end



numLayers=length(dh);
nx=zeros(1,numLayers);
ny=zeros(1,numLayers);
nz=zeros(1,numLayers);

layerStart=zeros(1,numLayers);
layerEnd=zeros(1,numLayers);

nGridpoints=0;
nGridpointsPerLayer=zeros(1,numLayers);


for ii=1:numLayers
    
    if  ii==numLayers || dh(ii)<dh(ii+1)
        layerEnd(ii)=interface(ii);
        ny(ii)=(layerEnd(ii)-layerStart(ii))/dh(ii)+1;
        if ii ~= numLayers
        layerStart(ii+1)=interface(ii)+dh(ii+1);
        end
    else
         layerEnd(ii)=interface(ii)-dh(ii);
         ny(ii)=(layerEnd(ii)-layerStart(ii))/dh(ii)+1;
         layerStart(ii+1)=interface(ii);
    end  
    
    nx(ii)=floor(NXmax/dh(ii));
    nz(ii)=floor(NZmax/dh(ii));
    if dh(ii)>1
        nx(ii)= nx(ii)+1;
        nz(ii)= nz(ii)+1;
    end

    nGridpointsPerLayer(ii)=nx(ii)*ny(ii)*nz(ii);
    
    nGridpoints=nGridpoints+nGridpointsPerLayer(ii);
    
end

% nx
% ny
% nz


values=zeros(1,nGridpoints);
size=zeros(1,nGridpoints);

for ii=1:nGridpoints
 %   ii=4998;
    index=ii-1;
    layer=0;
    %find out which in which layer the index is and reduce index to index
    %inside this layer
    for jj=1:numLayers
        index=index-nGridpointsPerLayer(jj);
        if index<0
            index=index+nGridpointsPerLayer(jj);
            layer=jj;
            break
        end
    end
    % calculate coordinates inside a subgrid
    

    coordinates.y(ii) = floor(index / (nx(layer) * nz(layer)));
    index = index - coordinates.y(ii) * (nx(layer) * nz(layer));
    
    coordinates.z(ii) = floor(index / (nx(layer)));
    index = index - coordinates.z(ii) * (nx(layer));
    
    coordinates.x(ii) = index;
    
    % coordinates in reference to the fine grid
    coordinates.x(ii)=coordinates.x(ii)*dh(layer);
    coordinates.y(ii)=coordinates.y(ii)*dh(layer);
    coordinates.z(ii)=coordinates.z(ii)*dh(layer); 
    
                
    % move the subgrid coordinates to global coordinates
    coordinates.y(ii)=coordinates.y(ii)+layerStart(layer);
   
    x=coordinates.x(ii);
    y=coordinates.y(ii);
    z=coordinates.z(ii);
    
    
    if (ii)==32441
        x
        y
        z
    end
    if (ii)==32237
        x
        y
        z
    end
    
%     if (ii)==3977
%         x
%         y
%         z
%     end
    
    
    % set values nd gridpointsize for plotting
    values(ii)=layer;
    size(ii)=dh(layer);
    
    
    %%test coordinate2index
 
    %find correct layer and set y coordinate to the y coordinate in the
    % subgrid to get local coordinates


    for jj=1:numLayers
        if coordinates.y(ii) <= layerEnd(jj) && coordinates.y(ii) >= layerStart(jj)
           testlayer=jj;          
           subGridCoordY=coordinates.y(ii)-layerStart(jj);
        end
    end

    % find index in the current subgrid
    testIndex= (((coordinates.x(ii))/dh(testlayer)) + (coordinates.z(ii)/dh(testlayer))*nx(layer) + (subGridCoordY/dh(testlayer))*nx(layer) * nz(layer));
    
    % add indices of the grids above the current subgrid to get the global
    % index
        for jj=2:testlayer
            testIndex=testIndex+nGridpointsPerLayer(jj-1);
        end          

    

    %$ test if coordinate2index worked
   % if ii ==2838
    if ii ~=testIndex+1
%          testlayer
%          subGridCoordY
%          coordinates.y(ii)
%           ii
%           testIndex
       break
    end
    
end

%% Plotting
% plot3(coordinates.x,coordinates.y,coordinates.z,'v')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% view(180,90)
%view(270,180)



scatter3(coordinates.x,coordinates.y,coordinates.z,size*50,values,'filled','s')
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)
%view(270,180)

% [xq,yq,zq] = meshgrid(0:1:NXmax, 0:1:100, 0:1:NZmax);
% vq = griddata(coordinates.x,coordinates.y,coordinates.z,values,xq,yq,zq);

% [xq,yq] = meshgrid(0:1:NXmax, 0:1:67);
% vq = griddata(coordinates.x,coordinates.y,values,xq,yq);

%
% figure(1)
% imagesc(vq)
% %mesh(vq)
% axis ij
