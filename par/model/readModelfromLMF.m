function [model]=readModelFromLMF(filename,NX,NY,NZ)
% readModelfromLMF  read 1D models into 3D array

%   model = readModelfromLMF(filename,NX,NY,NZ)
%   filename: filename of the model
%   NX,NY,NZ number of gridpoints in each direction+
%
%   major order of input should by X,Z,Y
%   output is model(Y,X,Z)

fileID = fopen(filename,'r');

% read the header [ID, indexType, valueType]
HEADER = fread(fileID, [1 3], 'int');
if HEADER(1) ~= hex2dec( '4711E01' )
   error('HEADER=%s not 4711E01 (dense vector)', dec2hex(HEADER(1)))
end
if HEADER(2) ~= 0
   error('data type for index values not int')
end
if HEADER(3) ~= 2
   error('data type for values not float')
end

NDIMS  = fread(fileID, [1 1], 'int');
if NDIMS~=1
    error('readModelFromLMF, NDIMS=%d must be 1', NDIMS);
end
SIZE = fread(fileID, [1 1], 'int');
if SIZE ~= NX * NY * NZ
    error('readModelFromLMF, SIZE=%d does not match NX=%d x NY=%d x NZ=%d', SIZE, NX, NY, NZ);
end

% now read the binary 'float' vector
A = fread(fileID, [1 SIZE], 'float');

% major order of vector is X, Z, Y
model=permute((reshape(A(:),[NX, NZ, NY])),[3 1 2]);

fclose(fileID);

end
