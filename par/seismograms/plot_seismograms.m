clearvars; close all;

filename='seismogram.mtx';

seismogram=read_seismogram(filename);

figure
plot(seismogram')
