function [V] = readVoltif(fname)
%Read tif image stack (3D)
%
%   [V] = readVoltif(fname)
%   -----------------------
%
%   Inputs:
%       > fname: file name
%   
%   Outputs:
%       >     V: image matrix
%
%   Keerthi Krishna PARVATHANENI  2018.02.01
%

info = imfinfo(fname);
nsl = numel(info);

A = imread(fname, 1);
V = zeros([size(A),nsl],class(A));

for i = 1:nsl
    V(:,:,i) = imread(fname,i);
end

end


