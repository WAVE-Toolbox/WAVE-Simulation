close all;
clear all
clc



%%
% Gridpoints of layer with dh=1 (every gridpoint)
NXmax=104;
NZmax=104;
NYmax=71;


interface(1)=4;
interface(2)=37;
interface(3)=41;
interface(4)=74;
interface(5)=76;

dh(1)=1;
dh(2)=3;
dh(3)=1;
dh(4)=3;
dh(5)=1;

numLayers=length(dh);

ny=zeros(1,length(dh));

lower=0;
higher=0;

for ii=1:numLayers
    
    if  ii==numLayers || dh(ii)<dh(ii+1)
        higher=interface(ii);
        ny(ii)=(higher-lower)/dh(ii)+1;
        if ii ~= numLayers
        lower=interface(ii)+dh(ii+1);
        end
    else
         higher=interface(ii)-dh(ii);
         ny(ii)=(higher-lower)/dh(ii)+1;
         lower=interface(ii);
    end  
    
end


ny
%%



DHmax=max(dh);
if (NXmax ~= floor(NXmax/DHmax)*DHmax+1 + floor(DHmax/2))
    message ="NXmax must be " + num2str(floor(NXmax/DHmax)*DHmax+1 + floor(DHmax/2));
    disp(message)
end


for ii=1:numLayers
    nx(ii)=floor(NXmax/dh(ii));
    nz(ii)=floor(NZmax/dh(ii));
    if dh(ii)>1
        nx(ii)= nx(ii)+1;
        nz(ii)= nz(ii)+1;
    end
end

%nGridpointsPerLayer=zeros(1,numLayers);
nGridpoints=0;
nGridpointsPerLayer=zeros(1,numLayers);

for ii=1:numLayers
    nGridpointsPerLayer(ii)=nx(ii)*ny(ii)*nz(ii);
    
    nGridpoints=nGridpoints+nGridpointsPerLayer(ii);
    
end

values=zeros(1,nGridpoints);
size=zeros(1,nGridpoints);

for ii=1:nGridpoints
    %ii=2148;
    layerindex=ii;
    layer=0;
    %find out which in which layer the index is
    for jj=1:numLayers
        layerindex=layerindex-nGridpointsPerLayer(jj);
        if layerindex<=0
            layerindex=layerindex+nGridpointsPerLayer(jj);
            layer=jj;
            break
        end
    end
    % calculate coordinates inside a subgrid
    
    values(ii)=layer;
    size(ii)=dh(layer);
    
    index=layerindex-1;
    coordinates.y(ii) = floor(index / (nx(layer) * nz(layer)));
    index = index - coordinates.y(ii) * (nx(layer) * nz(layer));
    
    coordinates.z(ii) = floor(index / (nx(layer)));
    index = index - coordinates.z(ii) * (nx(layer));
    
    coordinates.x(ii) = index;
    
    coordinates.x(ii)=coordinates.x(ii)*dh(layer);
    coordinates.y(ii)=coordinates.y(ii)*dh(layer);
    coordinates.z(ii)=coordinates.z(ii)*dh(layer);
    
    %     layer
    %     coordinates.y(ii)
    
    % move the subgrid coordinates to global coordinates
    if layer>1
        for jj=2:layer
            offset=ny(jj-1)*dh(jj-1)-dh(jj-1)+dh(jj) ;
            
            if dh(jj-1)>dh(jj)
                offset=offset+dh(jj-1)-dh(jj);
            end
            
            coordinates.y(ii)=coordinates.y(ii)+offset;
        end
    end
    
    
    
    %       coordinates.y(ii)
    
    %%test coordinate2index
    testlayer=1;
    subGridCoordY=coordinates.y(ii);
    
    %find correct layer and set y coordinate to the y coordinate in the
    % subgrid to get local coordinates

    for jj=2:numLayers
        nydummy=ny(jj-1)*dh(jj-1);
        if dh(jj-1) > 1
            nydummy=nydummy-dh(jj-1);
        end
        
        if subGridCoordY>nydummy
            
            testlayer=testlayer+1;
            
            subGridCoordY=subGridCoordY-ny(jj-1)*dh(jj-1)-dh(jj)+dh(jj-1);
            if dh(testlayer) < dh(testlayer-1)
               subGridCoordY=subGridCoordY - dh(jj-1)+dh(jj);
            end
   
        else
            break
        end
        
    end

    % find index in the current subgrid
    testIndex= (((coordinates.x(ii))/dh(testlayer)) + (coordinates.z(ii)/dh(testlayer))*nx(layer) + (subGridCoordY/dh(testlayer))*nx(layer) * nz(layer));
    
    % add indices of the grids above the current subgrid to get the global
    % index
    if testlayer > 1
        for jj=2:testlayer
            testIndex=testIndex+nGridpointsPerLayer(jj-1);
        end          
    end
    

    % test if coordinate2index worked
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
