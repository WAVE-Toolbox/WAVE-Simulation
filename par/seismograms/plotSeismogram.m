clearvars; close all;

addpath('../configuration')
config=conf('../configuration/configuration.txt');

%% Read seismogram
component='p';
shotnr=1;

format=config.getValue('FileFormat');
filename_base=config.getString('SeismogramFilename');
filename=['../',filename_base,'.shot_',num2str(shotnr),'.',component];
% readSeismogram 1=mtx 2=lmf
seismogram=readSeismogram(filename,format);

DT=config.getValue('seismoDT');

T=1*DT:DT:size(seismogram,2)*DT;

%% Plot seismogram
figure
for(trace=1:size(seismogram,1))
plot(T,seismogram(trace,:)/max(abs(seismogram(trace,:)))+trace,'black');
hold on
end
title('Normalized traces')
xlabel('Time in seconds')
ylabel('Traces')
axis([0 size(seismogram,2)*DT 0 size(seismogram,1)+1])