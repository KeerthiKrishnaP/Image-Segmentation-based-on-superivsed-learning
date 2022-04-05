function [K] = gaussianKernal(r,sig,ndim)
%compute the kernal matrix for window function
%
%   [K] = gaussianKernal(r,sig,ndim)
%   --------------------------------
%
%   Inputs:
%       >    r : kernal radius
%       >  sig : std of Gaussian distribution
%       > ndim : dimension (2 or 3)
%
%   Outputs:
%       >    K : Gaussian kernal
%
% Keerthi Krishna PARVATHANENI  2018.03.06
%

if ndim==2
    [y,x] = ndgrid(-r:1:r,-r:1:r);
    K = exp(-(x.^2+y.^2)./(2*sig^2)); %Gaussian average filter kernal
    K = K./(sum(K(:)));
    %figure;surf(x(:,:),y(:,:),w(:,:));set(gca,'ZLim',[0 max(w(:))])
elseif ndim==3
    [y,x,z] = ndgrid(-r:1:r,-r:1:r,-r:1:r);
    K = exp(-(x.^2+y.^2+z.^2)./(2*sig^2)); %Gaussian average filter kernal
    K = K./(sum(K(:)));
    %figure;surf(x(:,:,r+1),y(:,:,r+1),w(:,:,r+1));set(gca,'ZLim',[0 max(w(:))])
end
