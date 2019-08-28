clear all
close all
clc

%t=10^-3:10^-3:5;
out=load('seismogram.3D.acoustic.ci.p.mtx' ,'-ascii');

figure(1)
plot(out)

ref=load('seismogram.3D.acoustic.ref.p.mtx' ,'-ascii');
hold on

plot(-ref,'r')