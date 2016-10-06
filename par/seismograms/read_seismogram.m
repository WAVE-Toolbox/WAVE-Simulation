
function [Seismogram]=read_seismogram(filename)

fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE)

A=fscanf(fileID,'%d %d %e',[3 size(3)]);
Seismogram=accumarray([A(1,:)' A(2,:)'],A(3,:));

end