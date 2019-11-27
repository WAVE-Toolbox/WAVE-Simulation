function [model1D]=write3DModel2LMF(filename,model)
% write3DModel2LMF  writes 3D array models to 1D LMF files and returns
% model vector

%   write3DModel2LMF(filename,model)
%   filename: filename of the model
%   model: 3D array
%   input must be model(Y,X,Z)
%   major order of output is X,Z,Y
%   with X: fast dimension


model1D=reshape(permute(model,[2 3 1]),[],1);
writeVector2LMF(filename,model1D);

end
