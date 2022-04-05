function [divV] = div_5pts(varargin)
%compute the divergence of vectors whose components are defined in two 2D
%images or three 3D images.
% using 5-points numerical differentiation
%
%   [divV] = div_5pts(Vx,Vy)
%   [divV] = div_5pts(Vx,Vy,Vz)
%   ---------------------------
%
%   Inputs:
%       >    Vx,Vy : 2D images defining the x-,y-components of the vector
%       > Vx,Vy,Vz : 3D images defining the x-,y-,z-component of the vector
%
%   Outputs:
%       > divV : divergence image
%
% Keerthi Krishna PARVATHANENI  2018.02.27
%

narginchk(2,3);

if nargin==2
    [Vxx,Vxy] = imgrad_5pts(varargin{1});
    [Vyx,Vyy] = imgrad_5pts(varargin{2});
    divV = Vxx + Vyy;
elseif nargin==3
    [Vxx,Vxy,Vxz] = imgrad_5pts(varargin{1});
    [Vyx,Vyy,Vyz] = imgrad_5pts(varargin{2});
    [Vzx,Vzy,Vzz] = imgrad_5pts(varargin{3});
    divV = Vxx + Vyy + Vzz;
end