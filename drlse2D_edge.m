function [phi,Tdr,Ta,Te] = drlse2D_edge(phi,g,gx,gy,lambda,mu,alpha,eps,dt,varargin)
%compute the zero level set in a 2D image using DRLSE proposed by ...
%[Li et al., 2010]
%
%


% Notes:
%   1 - when calculating numerical derivation (grad, div, laplace),
%   numerical error is systematically added to the result, so if possible
%   avoid to perform div(grad) by using two numerical derivations. It is
%   strongly recommanded to use laplace instead of div(grad), so that only
%   one numerical derivation is employed!
%   Otherwise, the numerical results between del2(A) and div(grad(A)) could
%   be very different ! But why ? For now, it is not yet clear, but it
%   should be related to the numerical error in discrete derivations...
%
%   2 - the laplace operator is done by the MatLab function del2.m
%       Attention to the use of this function: L(A) = 2*N*del2(A), with N
%       the dimension of the input A.
%

narginchk(9,10);

% choose the potential function for distance regularization
if nargin==9
    pFct = 'double-well';
elseif nargin==10
    pFct = varargin{1};
end

% distance regularization term
[phix,phiy] = gradient(phi);
s = sqrt(phix.^2+phiy.^2);
s(s==0) = 1e-10;
phixs = phix ./ s;
phiys = phiy ./ s;
curvature=div(phixs,phiys);
if strcmp(pFct,'single-well')
    Tdr = 4*del2(phi)-curvature; 
elseif strcmp(pFct,'double-well')
    a = s>=0 & s<=1;
    b = s>1;
    ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);
    dps=((ps~=0).*ps+(ps==0)) ./ ((s~=0).*s+(s==0));
    Tdr = div(dps.*phix-phix,dps.*phiy-phiy) + 4*del2(phi);
end

% edge term
diracR_phi = diracR(phi,eps);
Te = diracR_phi .* ( (gx.*phixs+gy.*phiys) + g.*(curvature) );

% area term
Ta = g .* diracR_phi;

% update the level set function
phi = phi + dt * (mu*Tdr + lambda*Te + alpha*Ta);




%% internal functions
function [f] = diracR(x,eps)
%compute a regularized dirac function of a list of values
%
%   f = diracR(x,eps)
%   -----------------
%
%   Inputs:
%       >   x : (nxm matrix) list of values
%       > eps : regularization parameter
%
%   Outputs:
%       >   f : (nxm matrix) dirac function values
%
% Yang Chen 2.18.02.27
%

f = x.*0;
id = abs(x)<=eps;

eps2i = 1/2/eps;
pi_eps = pi/eps;
f(id) = eps2i .* (1 + cos(pi_eps.*x(id)));

