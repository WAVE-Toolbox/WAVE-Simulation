function [model]=readModelfromMtx(filename,NX,NY,NZ)
% readModelfromMtx  read 1D models into 3D array

%   model = readModelfromMtx(filename,NX,NY,NZ)
%   filename: filename of the model
%   NX,NY,NZ number of gridpoints in each direction+
%
%   major order of input should by X,Z,Y
%   output is model(Y,X,Z)



fileID = fopen(filename,'r');
HEADER = fgets(fileID);
SIZE = fgets(fileID);
size=str2num(SIZE);
A=fscanf(fileID,'%e',[1 size(1)]);
model=permute((reshape(A(:),[NX, NZ, NY])),[3 1 2]);
fclose(fileID);
end