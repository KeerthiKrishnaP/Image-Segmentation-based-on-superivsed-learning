function [S] = structTensor(A,r)
%Compute the structure tensor of an image
%    [S] = structTensor(A,r)
%    -----------------------
%    Inputs:
%       > A: 3-D image to be processed
%       > r: scalar, window "radius", the cubic window has a size of (2*r+1)^3
%
%    Output:
%       > S: structure data containing the six components of the structure
%            tensor and the window radius used
%
%    Keerthi Krishna PARVATHANENI  2018.01.25
%

% 5-point numerical differentiation ---------------------------------------
deriv5pt = [1 -8 0 8 -1];

Ix = convn(A,deriv5pt,'same');
Iy = convn(A,deriv5pt','same');
Iz = convn(A,permute(deriv5pt,[3,1,2]),'same');
% note: here we didnot yet carefully treat the image edges !
%       the convn uses zero-padded edges.
Ix = single(Ix);
Iy = single(Iy);
Iz = single(Iz);

% compute the local structure tensor --------------------------------------
S11 = Ix.^2;
S22 = Iy.^2;
S33 = Iz.^2;
S12 = Ix.*Iy;
S13 = Ix.*Iz;
S23 = Iy.*Iz;

% compute the average structure tensor ------------------------------------
% note: different operators are possible for this step.
%       - [Straumit2015] uses a uniform average filter;
%       - Here, we can also test with a Gaussian-type average filter (inactivated).
w = ones(2*r+1,2*r+1,2*r+1); %uniform average filter
    % sig = 1;
    % [y,x,z] = ndgrid(-r:1:r,-r:1:r,-r:1:r);
    % w = exp(-(x.^2+y.^2+z.^2)./(2*sig^2)); %Gaussian average filter
    %     figure;surf(x(:,:,r+1),y(:,:,r+1),w(:,:,r+1));set(gca,'ZLim',[0 1])
w = w./(sum(w(:)));
S11 = convn(S11,w,'same');
S22 = convn(S22,w,'same');
S33 = convn(S33,w,'same');
S12 = convn(S12,w,'same');
S13 = convn(S13,w,'same');
S23 = convn(S23,w,'same');

S = struct('S11',S11,'S22',S22,'S33',S33,'S12',S12,'S13',S13,'S23',S23,...
           'r',r);
