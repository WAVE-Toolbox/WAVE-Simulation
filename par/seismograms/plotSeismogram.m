clearvars; close all;

%% Read seismogram
filename='seismogram.shot_0.p';
% readSeismogram 1=mtx 2=lmf
seismogram=readSeismogram(filename,2);


DT=2e-3;


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