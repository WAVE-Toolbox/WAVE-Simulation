clearvars; close all;

%% Read seismogram
filename='seismogram.p.mtx';
seismogram=readSeismogram(filename);

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