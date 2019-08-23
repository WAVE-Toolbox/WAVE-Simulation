function writeModelToBinary(model,filename)

fileID = fopen(filename,'w');
fwrite(fileID,model,'float32','ieee-le');
fclose(fileID);
end

% Transform to su with
%suaddhead ns=NY < filename.bin | sushw key=gy,gx,dt a=0,0,DH b=DH,0,0 c=0,DH,0 j=NX,NX,0  > filename.su
