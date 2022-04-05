function [A] = addMarg(A,dim,val,len)
%Add margin layers aournd a 3D image
%
%   [A] = addMarg(A,dim,val)
%
%   Inputs:
%       >   A : input 3D image
%       > dim : vector(length<=3) defining the dimensions along which the
%               margin layers are added 
%       > val : vector(length<=3) values to be given to the added margin
%               layers along every dimension
%       > len : vector(length<=3) lengths of the margins to be added
%
%   Outputs:
%       >   A : output 3D image
%
% Keerthi Krishna PARVATHANENI 2018.05.17
% modified Keerthi Krishna PARVATHANENI 2018.10.15 (only one layer in each dimension)
%

if ismember(1,dim)
    id = dim==1;
    sc = val(id);
    siz = size(A);
    tmp = permute(zeros([siz([2,3]),len(id)]),[3 1 2]) + sc;
    %A = cat(1,tmp,A,tmp);
    A = cat(1,tmp,A);
end

if ismember(2,dim)
    id = dim==2;
    sc = val(id);
    siz = size(A);
    tmp = permute(zeros([siz([1,3]),len(id)]),[1 3 2]) + sc;
    %A = cat(2,tmp,A,tmp);
    A = cat(2,tmp,A);
end

if ismember(3,dim)
    id = dim==3;
    sc = val(id);
    siz = size(A);
    tmp = permute(zeros([siz([1,2]),len(id)]),[1 2 3]) + sc;
    %A = cat(3,tmp,A,tmp);
    A = cat(3,tmp,A);
end


