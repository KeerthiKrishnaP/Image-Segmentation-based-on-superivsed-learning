function [varargout] = imgrad_5pts(A)
%compute the gradiant of a 2D/3D image
% using 5-point numerical differentiation
%
%	[Ax,Ay] = imgrad_5pts(A)
%	[Ax,Ay,Az] = imgrad_5pts(A)
%   ---------------------------
%
%   Inputs:
%       > A : 2D / 3D image to be processed
%
%   Outputs:
%       > Ax,Ay,Az: gradiant components in x,y,z directions
%
% note: here we didnot yet carefully treat the image edges !
%       the convn uses zero-padded edges.
%
% Keerthi Krishna PARVATHANENI  2018.02.27
%


deriv5pt = [1 -8 0 8 -1];

siz = size(A);

if length(siz)==2
    Ax = convn(A,deriv5pt,'same');
    Ay = convn(A,deriv5pt','same');
    varargout{1} = Ax;
    varargout{2} = Ay;
elseif length(siz)==3
    Ax = convn(A,deriv5pt,'same');
    Ay = convn(A,deriv5pt','same');
    Az = convn(A,permute(deriv5pt,[3,1,2]),'same');
    varargout{1} = Ax;
    varargout{2} = Ay;
    varargout{3} = Az;
end
