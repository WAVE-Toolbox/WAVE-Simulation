
function [Seismogram]=readSeismogram(filename,fileFormat)
% readSeismogram  read 1D models into 3D array

%   Seismogram = readSeismogram(filename,fileFormat)
%   filename: filename of the seismogram without suffix
%   fileFormat: 1=mtx 2=lmf
%
%   output is Seismogram(ntr,ns);

%mtx
if fileFormat==1
    if ~contains(filename,'.mtx')
        filename=[filename '.mtx'];
    end
    fileID = fopen(filename,'r');
    HEADER = fgets(fileID);
    SIZE = fgets(fileID);
    SIZE=str2num(SIZE);
    Seismogram=fscanf(fileID,'%e',[ SIZE(1) SIZE(2)]);
end

%lmf
if fileFormat==2

    % read the header [ID, indexType, valueType]
    if ~contains(filename,'.lmf')
        filename=[filename '.lmf']; 
    end
    fileID = fopen(filename,'r');

    %read the header [ID, indexType, valueType]
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
    if NDIMS~=2
        error('readModelFromLMF, NDIMS=%d must be 2', NDIMS);
    end

    SIZE = fread(fileID, [1 2], 'int');

    A = fread(fileID, [SIZE(2) SIZE(1)], 'float');
    fclose(fileID);

    Seismogram=A';    
end

