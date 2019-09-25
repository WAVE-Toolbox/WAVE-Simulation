function writeVector2LMF(filename,vector)

%filename
fileID = fopen(filename,'w');

%HEADER
HEADER=[hex2dec( '4711E01' ) 0 2];
fwrite(fileID,HEADER,'int');

%NDIMS
NDIMS=1;
fwrite(fileID,NDIMS,'int');

%SIZE
SIZE=length(vector);
fwrite(fileID,SIZE,'int');

fwrite(fileID,vector,'float32','ieee-le');
fclose(fileID);


end