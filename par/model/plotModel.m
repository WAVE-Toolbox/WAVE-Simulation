clearvars; close all;

addpath('../configuration')
config=conf('../configuration/configuration.txt');

%% Define input parameter
comp='vp';

format=config.getValue('FileFormat');
if format==1
    format_str='mtx';
elseif format==2
    format_str='lmf';
end
filename_base=config.getString('ModelFilename'); % File name of the model
filename=['../',filename_base,'.',comp,'.',format_str];
NX=config.getValue('NX');  % Number of grid points in X
NY=config.getValue('NY');  % Number of grid points in Y
NZ=config.getValue('NZ');  % Number of grid points in Z
DH=config.getValue('DH');   % Spatial grid sampling
LAYER=1; % Define layer of 3D model to display as 2D slice

%% Read model
if strcmp(format_str,'mtx')
    model=readModelfromMtx(filename,NX,NY,NZ);
elseif strcmp(format_str,'lmf')
    model=readModelfromLMF(filename,NX,NY,NZ);
else
    disp('unknown format to read model')
end
X=0:DH:(NX*DH-DH);
Y=0:DH:(NY*DH-DH);

%% Plot
figure
imagesc(X,Y,squeeze(model(:,:,LAYER)))
colorbar
xlabel('X in meter')
ylabel('Y in meter')