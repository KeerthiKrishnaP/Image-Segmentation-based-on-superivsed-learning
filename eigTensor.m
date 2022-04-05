function [varargout] = eigTensor(S)
%Compute the eigenthings of the structure tensor
%
%    [lambda,vec,beta] = eigTensor(S)
%    [beta] = eigTensor(S)
%    --------------------------------
%    Inputs:
%       >      S: structure data containing the six components of the
%                 structure tensor.
%
%    Outputs:
%       > lambda: eigenvalues of the structure tensor (lambda1<lambda2<lambda3)
%       >    vec1: eigenvector corresponding to lambda1
%       >    vec2: eigenvector corresponding to lambda2
%       >    vec3: eigenvector corresponding to lambda3
%       >   beta: anisotropy parameter
%
%
% Keerthi Krishna PARVATHANENI 2018.01.25
% Keerthi Krishna PARVATHANENI 2018.03.02 (version chapo* using "mex eig3volume_YC.c")
% eig3volume_YC-- MEX file to compute Structre tensor developed along with YANG Chen 

% eigenthings
[L1,L2,L3,V1x,V1y,V1z,V2x,V2y,V2z,V3x,V3y,V3z] = ...
    eig3volume_YC(S.S11(:),S.S12(:),S.S13(:),S.S22(:),S.S23(:),S.S33(:));

lambda = [L1(:),L2(:),L3(:)]; clear L1 L2 L3
vec1 = [V1x(:),V1y(:),V1z(:)]; clear V1x V1y V1z
vec2 = [V2x(:),V2y(:),V2z(:)]; clear V2x V2y V2z
vec3 = [V3x(:),V3y(:),V3z(:)]; clear V3x V3y V3z

% nvx = length(S.S11(:));
% lambda = zeros(nvx,3,'single');
% vec1 = zeros(nvx,3,'single');
% vec2 = zeros(nvx,3,'single');
% vec3 = zeros(nvx,3,'single');
% for ivx = 1:nvx
%     S_ivx = [S.S11(ivx),S.S12(ivx),S.S13(ivx);
%              S.S12(ivx),S.S22(ivx),S.S23(ivx);
%              S.S13(ivx),S.S23(ivx),S.S33(ivx)];
%     [vec_i,val_i] = eig(S_ivx);
%     lambda(ivx,:) = [val_i(1,1) val_i(2,2) val_i(3,3)];
%     vec1(ivx,:) = vec_i(:,1);
%     vec2(ivx,:) = vec_i(:,2);
%     vec3(ivx,:) = vec_i(:,3);
% end

% anisotropy parameter: beta
beta = 1 - lambda(:,1)./lambda(:,3);
beta(lambda(:,3)==0) = 0;

% outputs
nargoutchk(1,5);
if nargout==1
    varargout{1} = beta;
elseif nargout==2
    varargout{1} = lambda;
    varargout{2} = vec1;
elseif nargout==3
    varargout{1} = lambda;
    varargout{2} = vec1;
    varargout{3} = beta;
elseif nargout==5
    varargout{1} = lambda;
    varargout{2} = vec1;
    varargout{3} = vec2;
    varargout{4} = vec3;
    varargout{5} = beta;
end

