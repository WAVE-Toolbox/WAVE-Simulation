close all;
clearvars
clc





%% Input
% Gridpoints of layer with dh=1 (every gridpoint)
NXmax=305;
NZmax=1;
NYmax=305;

interface(1)=100;
interface(2)=151;
interface(3)=160;


dh(1)=1;
dh(2)=3;
dh(3)=1;
dh(4)=3;



interface(length(interface)+1)=NYmax-1;
numLayers=length(dh);

%%
DHmax=max(dh);
while true
  if (NXmax == floor(NXmax/DHmax)*DHmax+1 + floor(DHmax/2))   
      break
  end
  NXmax=NXmax-1;
 disp(['NX varGrid has been changed from ' num2str(NXmax+1) ' to ' num2str(NXmax)]);
end

if NZmax~=1
while true
  if (NZmax == floor(NZmax/DHmax)*DHmax+1 + floor(DHmax/2))   
      break
  end
  NZmax=NZmax-1;
 disp(['NZ varGrid has been changed from ' num2str(NZmax+1) ' to ' num2str(NZmax)]);
end
end

ii=1;
%check interfaces
while ii<numLayers
    %test
    ii=ii+1;
    test=(interface(ii)-interface(ii-1))/dh(ii);
    
    if floor(test)~=test     
        interface(ii)=interface(ii)-1;
        disp(['interface ' num2str(ii) ' has been moved from ' num2str(interface(ii)+1) ' to ' num2str(interface(ii))]);
        ii=ii-1;
    end
    
    
end
%%


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


values=zeros(1,nGridpoints);
size=zeros(1,nGridpoints);

for ii=1:nGridpoints
    index=ii-1;

    
    [coordinates.x(ii),coordinates.y(ii),coordinates.z(ii),layer]=index2coordinate(index,nGridpointsPerLayer,layerStart,dh,nx,nz); 
    
    % set values nd gridpointsize for plotting
    values(ii)=layer;
    size(ii)=dh(layer);
    
    %%test coordinate2index
 
    [testIndex, testLayer] = coord2index(coordinates.x(ii),coordinates.y(ii),coordinates.z(ii),nGridpointsPerLayer,layerStart,layerEnd,dh,nx,nz);
    
    if layer~=testLayer || index~=testIndex
        display(testIndex);
        display(index);
        error('something wrong')
    end
     
    
end




%% Plotting

scatter3(coordinates.x,coordinates.y,coordinates.z,size*50,values,'filled','s')
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)
view(0,90)

% [xq,yq,zq] = meshgrid(0:1:NXmax, 0:1:100, 0:1:NZmax);
% vq = griddata(coordinates.x,coordinates.y,coordinates.z,values,xq,yq,zq);

% [xq,yq] = meshgrid(0:1:NXmax, 0:1:67);
% vq = griddata(coordinates.x,coordinates.y,values,xq,yq);

%
% figure(1)
% imagesc(vq)
% %mesh(vq)
% axis ij
%%


function [x,y,z,layer] = index2coordinate(index,nGridpointsPerLayer,layerStart,dh,nx,nz)
    
    for jj=1:length(dh)
        index=index-nGridpointsPerLayer(jj);
        if index<0
            index=index+nGridpointsPerLayer(jj);
            layer=jj;
            break
        end
    end
    % calculate coordinates inside a subgrid
    

    y = floor(index / (nx(layer) * nz(layer)));
    index = index - y * (nx(layer) * nz(layer));
    
    z = floor(index / (nx(layer)));
    index = index - z * (nx(layer));
    
    x = index;
    
    % coordinates in reference to the fine grid
    x=x*dh(layer);
    y=y*dh(layer);
    z=z*dh(layer); 
    
                
    % move the subgrid coordinates to global coordinates
    y=y+layerStart(layer);
   
end


function [index,testlayer] = coord2index(x,y,z,nGridpointsPerLayer,layerStart,layerEnd,dh,nx,nz)

if x>=max(nx) 
    error('nx out of bounds')
end
if z>=max(nz) 
    error('nz out of bounds')
end
testlayer=1;
subGridCoordY=0;

    for jj=1:length(dh)
        if y <= layerEnd(jj) && y >= layerStart(jj)
           testlayer=jj;          
           subGridCoordY=y-layerStart(jj);
        end
    end

    % find index in the current subgrid
    index= floor(x/dh(testlayer)) + (z/dh(testlayer))*nx(testlayer) + (subGridCoordY/dh(testlayer))*nx(testlayer) * nz(testlayer);
    
    % add indices of the grids above the current subgrid to get the global
    % index
        for jj=2:testlayer
            index=index+nGridpointsPerLayer(jj-1);
        end  

end