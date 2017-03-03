clearvars; close all;

%% Read seismogram
filename='seismogram.p.mtx';
seismogram=readSeismogram(filename);

%% Plot seismogram
figure
for(trace=1:size(seismogram,1))
plot(seismogram(trace,:)/max(abs(seismogram(trace,:)))+trace,'black');
hold on
end
title('Normalized traces')
xlabel('Samples')
ylabel('Traces')
axis([0 size(seismogram,2) 0 size(seismogram,1)+1])